#!/usr/bin/env python3
"""Destructive Make checks confined to a newly created external scratch root."""
from pathlib import Path
import csv,hashlib,json,os,re,shutil,subprocess,sys,time,signal
repo=Path(__file__).resolve().parents[1];root=Path(sys.argv[1]).resolve()
assert not root.is_relative_to(repo)
root.mkdir() # Never reuse or clean a user's directory.
env={**os.environ,'GIT_OPTIONAL_LOCKS':'0','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
results=[]
def run(args,cwd,name,timeout=240):
 log=root/(name+'.log');start=time.monotonic()
 with log.open('w')as output:
  process=subprocess.Popen(args,cwd=cwd,env=env,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=process.wait(timeout=timeout)
  except subprocess.TimeoutExpired:
   os.killpg(process.pid,signal.SIGTERM);process.wait();code=124
 results.append([name,code,round(time.monotonic()-start,3),' '.join(map(str,args)),str(cwd)])
 with(root/'results.csv').open('w')as f:
  w=csv.writer(f);w.writerow(['check','exit','seconds','command','cwd']);w.writerows(results)
 print(name,code,flush=True)
 return code

def copy_source(dest):
 dest.mkdir()
 for name in ['Makefile.am','Makefile.in','configure','configure.ac','aclocal.m4','config.h.in','compile','install-sh','missing','depcomp','config.guess','config.sub']:
  p=repo/name
  if p.exists():shutil.copy2(p,dest/name)
 for directory in ['src','distribution','interpolation']:
  (dest/directory).mkdir()
  for p in (repo/directory).iterdir():
   if p.is_file() and (p.suffix=='.f90' or p.name in ['Makefile.am','Makefile.in','README']):shutil.copy2(p,dest/directory/p.name)
 return dest

def artifacts(p):return {str(q.relative_to(p)):q.stat().st_mtime_ns for q in p.rglob('*.o')}
def generated(p):return sorted(str(q.relative_to(p))for ext in ['*.o','*.mod','*.smod']for q in p.rglob(ext))
def no_refresh(name):
 text=(root/(name+'.log')).read_text()
 assert not re.search(r'am--refresh|aclocal-1|automake-1|config.status --recheck|&& autoconf',text),name

serial=copy_source(root/'serial')
assert run(['./configure'],serial,'serial_configure')==0
assert run(['make','-j1'],serial,'serial_build')==0
assert (serial/'src/ALPS').is_file()
parallel=copy_source(root/'parallel')
fixture=parallel/'unrelated_fixture';fixture.mkdir()
(fixture/'obsolete.f90').write_text('module alps_analyt\ninteger :: obsolete_interface\nend module\n')
assert run([env.get('FC','mpifort'),'-c','obsolete.f90'],fixture,'unrelated_stub')==0
(fixture/'protected.in').write_text('scientific fixture sentinel\n')
protected={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in fixture.iterdir() if p.is_file()}
assert run(['./configure'],parallel,'parallel_configure')==0
assert run(['make','-j4'],parallel,'parallel_build')==0
no_refresh('parallel_build')
# Clean twice, checking every generated artifact and protected fixture each time.
for n in [1,2]:
 for d in ['src','distribution']:(parallel/d/'owned_cleanup.smod').write_text('owned test submodule artifact\n')
 assert run(['make','clean'],parallel,'clean_'+str(n))==0
 assert not [p for p in generated(parallel) if not p.startswith('unrelated_fixture/')],generated(parallel)
 assert all(hashlib.sha256((fixture/name).read_bytes()).hexdigest()==value for name,value in protected.items())
 assert run(['make','-j4'],parallel,'rebuild_'+str(n))==0
 no_refresh('rebuild_'+str(n))

incremental=[]
for provider,expected in [
 ('src/ALPS_analyt.f90',{'src/ALPS_analyt.o','src/ALPS_com.o','src/ALPS_fns_rel.o','src/ALPS_fns.o','src/ALPS.o'}),
 ('src/ALPS_rel_cache.f90',{'src/ALPS_rel_cache.o','src/ALPS_fns_rel.o','src/ALPS_fns.o','src/ALPS.o'}),
 ('distribution/distribution_analyt.f90',{'src/distribution_analyt.o','src/ALPS_analyt.o','src/ALPS_com.o','src/ALPS_fns_rel.o','src/ALPS_fns.o','src/ALPS.o','distribution/distribution_analyt.o','distribution/generate_distribution.o'})]:
 before=artifacts(parallel);os.utime(parallel/provider,None)
 name='incremental_'+Path(provider).stem
 assert run(['make','-j4'],parallel,name)==0
 after=artifacts(parallel);changed={p for p in after if after[p]!=before[p]}
 assert changed==expected,(provider,changed,expected)
 no_refresh(name);incremental.append({'provider':provider,'rebuilt_objects':sorted(changed)})
(root/'incremental.json').write_text(json.dumps(incremental,indent=2)+'\n')
assert run(['make','-j4'],parallel,'noop')==0
assert not re.search(r' -c -o ',(root/'noop.log').read_text());no_refresh('noop')
# Build directory ownership: building distribution first must not deprive src of
# its module; a clean in either directory must not delete the other's artifacts.
assert run(['make','clean'],parallel,'ownership_clean')==0
assert run(['make','-C','distribution','-j4'],parallel,'distribution_first')==0
owned=parallel/'distribution/distribution_analyt.o';mtime=owned.stat().st_mtime_ns
assert run(['make','-C','src','-j4'],parallel,'src_after_distribution')==0
assert owned.stat().st_mtime_ns==mtime
assert run(['make','-C','src','clean'],parallel,'src_only_clean')==0
assert owned.is_file() and (parallel/'distribution/alps_distribution_analyt.mod').is_file()
assert run(['make','-j4'],parallel,'ownership_rebuild')==0
# Two unconfigured-source out-of-tree builds; stale source-directory interfaces
# are deliberately left in place. Only build-local interfaces may be consumed.
source=copy_source(root/'source')
shutil.copy2(fixture/'alps_analyt.mod',source/'src/alps_analyt.mod')
source_stale_hash=hashlib.sha256((source/'src/alps_analyt.mod').read_bytes()).hexdigest()
for name in ['out_a','out_b']:
 build=root/name;build.mkdir()
 assert run([str(source/'configure')],build,name+'_configure')==0
 assert run(['make','-j4'],build,name+'_build')==0
 assert run(['make','clean'],build,name+'_clean')==0
 assert not generated(build)
 assert hashlib.sha256((source/'src/alps_analyt.mod').read_bytes()).hexdigest()==source_stale_hash
 assert run(['make','-j4'],build,name+'_rebuild')==0
# A real maintainer edit requires regeneration once, not on every object target.
os.utime(parallel/'src/Makefile.am',None)
assert run(['make','-C','src','ALPS_com.o'],parallel,'am_edit_refresh')==0
assert run(['make','-C','src','ALPS_com.o'],parallel,'am_edit_noop')==0
no_refresh('am_edit_noop')
os.utime(parallel/'configure.ac',None)
assert run(['make','-C','src','ALPS_com.o'],parallel,'configure_edit_refresh')==0
assert run(['make','-C','src','ALPS_com.o'],parallel,'configure_edit_noop')==0
no_refresh('configure_edit_noop')
# Diagnostic only: -B forces Makefile regeneration as well as object compilation.
# Its effects stay in this disposable tree; it is never a normal build recipe.
run(['make','-C','src','-n','-B','ALPS_analyt.o'],parallel,'forced_refresh_probe',timeout=20)
(root/'source_hashes.json').write_text(json.dumps({str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest()for d in ['src','distribution','interpolation']for p in (repo/d).iterdir()if p.is_file()and(p.suffix in ['.am','.in','.f90'])},indent=2)+'\n')
print('PASS: clean serial/parallel, scoped cleanup, incremental graph, independent module ownership and out-of-tree builds',flush=True)

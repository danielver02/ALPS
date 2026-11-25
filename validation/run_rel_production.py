#!/usr/bin/env python3
"""Ordinary production activation checks; copied-source observation only."""
from pathlib import Path
import sys,shutil,subprocess,os,re,csv,json,hashlib,difflib
repo=Path(__file__).resolve().parents[1];root=Path(sys.argv[1]).resolve()
assert not root.is_relative_to(repo);root.mkdir();(root/'source').mkdir();(root/'build').mkdir();(root/'cache').mkdir()
for p in(repo/'src').glob('*.f90'):shutil.copy2(p,root/'source'/p.name)
for p in [repo/'distribution/distribution_analyt.f90',repo/'validation/make_rel_susceptibility_case.f90',Path(__file__).with_name('rel_production_diagnostics.f90')]:shutil.copy2(p,root/'source'/p.name)
p=root/'source/ALPS.f90';s=p.read_text()
for marker in ['     call initialize_rel_cheb(rel_cheb_status)','  call complete_nr_cheb_initialization(nr_cheb_status)','  call startup_report','  call mpi_finalize (ierror)']:
 assert s.count(marker)==1,('production-main observation marker changed',marker)
s=s.replace('     call initialize_rel_cheb(rel_cheb_status)','     call initialize_rel_cheb(rel_cheb_status)\n     call r1_dump',1).replace('  call complete_nr_cheb_initialization(nr_cheb_status)','  call complete_nr_cheb_initialization(nr_cheb_status)\n  call r1_stats(\"initial\")',1);s=s.replace('  call startup_report','  call startup_report\n  call r2_preflight',1);s=s.replace('  call mpi_finalize (ierror)','  call r1_stats(\"final\")\n  call mpi_finalize (ierror)',1);p.write_text(s)
p=root/'source/ALPS_fns.f90';s=p.read_text()
marker='    call mpi_bcast(disp, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierror)'
assert s.count(marker)==1;s=s.replace(marker,'    if(proc0) call r2_tensor(om,disp)\n'+marker)
marker='     tmpDisp=disp(wroots(iw))';assert s.count(marker)==1
s=s.replace(marker,marker+'\n     call r2_root_record(wroots(iw),tmpDisp,iflag)');p.write_text(s)
flags=['-O2','-g','-fcheck=bounds','-fbacktrace','-ffree-line-length-none','-ffpe-trap=invalid,zero,overflow']
names=['ALPS_rel_cache','ALPS_var','distribution_analyt','ALPS_io','ALPS_analyt','ALPS_com','ALPS_NHDS','ALPS_fns_rel','ALPS_fns','ALPS_nyquist','ALPS_nyquist_growth','ALPS_nyquist_driver','rel_production_diagnostics']
for n in names:subprocess.run(['mpifort',*flags,'-c',str(root/'source'/(n+'.f90'))],cwd=root/'build',check=True)
subprocess.run(['mpifort',*flags,'-o','ALPS',*[n+'.o'for n in names],str(root/'source/ALPS.f90'),'-llapack','-lblas'],cwd=root/'build',check=True)
subprocess.run(['mpifort',*flags,'-o','make_case',str(root/'source/make_rel_susceptibility_case.f90')],cwd=root/'build',check=True)
subprocess.run(['mpifort',*flags,'-o','ALPS_plain',*[n+'.o'for n in names],str(repo/'src/ALPS.f90'),'-llapack','-lblas'],cwd=root/'build',check=True)
subprocess.run(['mpifort',*flags,'-o','support',*[n+'.o'for n in names],str(repo/'validation/test_rel_source_support.f90'),'-llapack','-lblas'],cwd=root/'build',check=True)
with(root/'support.log').open('w') as f:subprocess.run([str(root/'build/support')],stdout=f,check=True)
p=root/'source/ALPS_com.f90';normal=p.read_text()
marker='    if(.not.proc0) call rel_cheb_import(header(2),header(3),metadata,mapping,coefficients,local_status)'
assert normal.count(marker)==1
p.write_text(normal.replace('use alps_var, only : proc0,ierror','use alps_var, only : proc0,ierror,iproc').replace(marker,'    if(iproc==1) mapping(2,1)=-1.d0 ! validation-only worker import corruption\n'+marker))
subprocess.run(['mpifort',*flags,'-c',str(p)],cwd=root/'build',check=True)
subprocess.run(['mpifort',*flags,'-o','ALPS_badimport',*[n+'.o'for n in names],str(root/'source/ALPS.f90'),'-llapack','-lblas'],cwd=root/'build',check=True)
p.write_text(normal)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
cases=['method2','writeout','cached_fresh','cached_hit','logfalse','mixed_bm','mixed_nr','mixed_rel_method1','no_selection','bm_relflag','source_zero','source_negative','source_nan','source_missing','source_short','invalid_degree','invalid_kind','worker_import_failure']
with(root/'results.csv').open('w')as f:
 w=csv.writer(f);w.writerow(['case','exit','map_exists','row_diagnostics','initialization_messages']);f.flush()
 for case in cases:
  p=root/case;(p/'distribution').mkdir(parents=True);(p/'solution').mkdir()
  subprocess.run([str(root/'build/make_case'),str(repo),'coarse','6','24'],cwd=p,check=True)
  s=(p/'case.in').read_text()
  s=re.sub(r'(?im)^AC_method=.*','AC_method=2',s);s=re.sub(r'(?im)^order=.*','order=2',s)
  for k,v in [('use_map','T'),('determine_minima','F'),('omi','.04d0'),('omf','.05d0'),('gami','-.002d0'),('gamf','-.001d0'),('ni','2'),('nr','2')]:s=re.sub(r'(?im)^'+k+r'\s*=.*',k+'='+v,s)
  if case in ['cached_fresh','cached_hit']:
   s=re.sub(r"(?im)^rel_cache=.*","rel_cache='auto'",s)
   s=re.sub(r"(?im)^rel_cache_dir=.*","rel_cache_dir='"+str(root/'cache')+"'",s)
  if case=='no_selection':s=re.sub(r'(?im)^AC_method=.*','AC_method=1',s)
  if case=='invalid_kind':s=s.replace('kind=1','kind=2',1)
  if case=='writeout':s=re.sub(r'(?im)^writeOut=.*','writeOut=T',s)
  if case=='logfalse':s=re.sub(r'(?im)^log_fit=.*','log_fit=F',s)
  if case=='invalid_degree':s=s.replace('order=2','order=-1',1)
  a=s.index('&spec_2');b=a+re.search(r'(?m)^/\s*$',s[a:]).start();part=s[a:b]
  if case in ['mixed_bm','bm_relflag']:part=re.sub(r'(?im)^use_bM=.*','use_bM=T',part)
  if case in ['mixed_bm','mixed_nr']:part=re.sub(r'(?im)^relat=.*','relat=F',part)
  if case=='mixed_rel_method1':part=re.sub(r'(?im)^AC_method=.*','AC_method=1',part)
  s=s[:a]+part+s[b:]
  if case=='mixed_rel_method1':
   a=s.index('&ffit_2_1');b=a+re.search(r'(?m)^/\s*$',s[a:]).start();part=s[a:b];part=re.sub(r'(?im)^fit_type_in=.*','fit_type_in=4',part);s=s[:a]+part+s[b:]
  (p/'case.in').write_text(s)
  if case in ['source_zero','source_negative','source_nan','source_short','source_missing']:
   q=p/'distribution/susceptibility.1.array';lines=q.read_text().splitlines();values=lines[0].split();values[-1]={'source_zero':'0.0','source_negative':'-1.0','source_nan':'NaN','source_short':'0.0','source_missing':'0.0'}[case];lines[0]=' '.join(values);q.write_text('\n'.join(lines)+'\n')
   if case=='source_missing':q.unlink()
   if case=='source_short':q.write_text('\n'.join(lines[:10])+'\n')
  with(p/'stdout.log').open('w')as out,(p/'stderr.log').open('w')as err:
   try:result=subprocess.run(['mpirun','--oversubscribe','-np','4',str(root/'build'/('ALPS_badimport' if case=='worker_import_failure' else 'ALPS')),'case.in'],cwd=p,env=env,stdout=out,stderr=err,timeout=120);code=result.returncode
   except subprocess.TimeoutExpired:code=124
  text=(p/'stdout.log').read_text();w.writerow([case,code,(p/'solution/case.map').exists(),(p/'r1_rows.csv').exists(),text.count('R1_STATS')]);f.flush();print(case,code,flush=True)
(root/'source_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(root/'source').glob('*.f90')},indent=2)+'\n')

# Same arithmetic against preserved pre-activation ordinary-main observations.
archive=repo/'audits/rel_cheb_activation_r1'
for name in ['r1_nodes.csv','r1_rows.csv','r1_bounds.csv']:
 assert (root/'method2'/name).read_bytes()==(archive/name).read_bytes(),name
for case in ['method2','mixed_bm','mixed_nr','mixed_rel_method1']:
 assert (root/case/'solution/case.map').read_bytes()==(archive/(case+'_map.txt')).read_bytes(),case
assert (root/'method2/solution/case.map').read_bytes()==(root/'writeout/solution/case.map').read_bytes()
assert (root/'cached_fresh/solution/case.map').read_bytes()==(root/'cached_hit/solution/case.map').read_bytes()
assert 'origin=cache' in (root/'cached_hit/stdout.log').read_text()
assert not list((root/'writeout/distribution').glob('*.poly_fit_result*'))
# Prove the plain, uninstrumented ordinary main produces the same map.
p=root/'plain';shutil.copytree(root/'method2/distribution',p/'distribution');(p/'solution').mkdir();shutil.copy2(root/'method2/case.in',p/'case.in')
with(p/'stdout.log').open('w') as f:
 subprocess.run(['mpirun','--oversubscribe','-np','4',str(root/'build/ALPS_plain'),'case.in'],cwd=p,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=120,check=True)
assert (p/'solution/case.map').read_bytes()==(root/'method2/solution/case.map').read_bytes()
rows=list(csv.DictReader((root/'results.csv').open()))
for row in rows:
 expected=0 if row['case'] in ['method2','writeout','cached_fresh','cached_hit','mixed_bm','mixed_nr','mixed_rel_method1','no_selection'] else None
 if expected==0:assert row['exit']=='0',row
 else:assert int(row['exit']) not in [0,124,136],row
print('PASS: archived models/maps, plain production, output, fresh/cache, mixed paths and explicit failures',flush=True)

# Independent Python geometry counts on the same two bounded fixture grids.
import math
support=(root/'support.log').read_text()
observed=re.findall(r'cone_nodes=\s*(\d+) inside=\s*(\d+) exterior=\s*(\d+)',support)
for index,(rmin,rmax,umin,umax) in enumerate([(0.,5.,-5.,5.),(.4,1.2,-.7,1.3)]):
 counts=[0,0,0];gmax=math.sqrt(1+rmax*rmax+max(umin*umin,umax*umax))
 for i in range(51):
  g=1+(gmax-1)*i/50
  for j in range(201):
   u=umin+(umax-umin)*j/200
   if g*g-1<u*u:continue
   r=math.sqrt(max(0.,g*g-1-u*u));tol=128*sys.float_info.epsilon*max(1.,rmax,abs(umin),abs(umax),r,abs(u))
   inside=rmin-tol<=r<=rmax+tol and umin-tol<=u<=umax+tol
   counts[0]+=1;counts[1]+=inside;counts[2]+=not inside
 assert tuple(counts)==tuple(map(int,observed[index])),(counts,observed[index])
assert observed[0]==('7217','4978','2239')
assert observed[2]==('4','4','0'),observed[2]
print('PASS: independent source bounds, lower-perpendicular/asymmetric domain and boundary tolerance',flush=True)

comparisons=[]
for kind in ['iso','aniso']:
 archived=repo/'audits/task3e_roots'/(kind+'_np4_method2')
 reference=next(csv.DictReader((archived/'roots.csv').open()))
 tensors=list(csv.DictReader((archived/'tensor.csv').open()))
 for output in [False,True]:
  p=root/(kind+'_root_'+str(output));(p/'distribution').mkdir(parents=True);(p/'solution').mkdir()
  subprocess.run([str(root/'build/make_case'),str(repo),kind,'50','200'],cwd=p,check=True)
  s=(p/'case.in').read_text()
  replacements={'AC_method':'2','order':'2','nroots':'1','use_map':'F','n_scan':'0','secant_method':'1','numiter':'150','D_tol':'1.d-11','writeOut':'T' if output else 'F'}
  for k,v in replacements.items():s=re.sub(r'(?im)^'+k+r'\s*=.*',k+'='+v,s)
  a=s.index('&guess_1');b=a+re.search(r'(?m)^/\s*$',s[a:]).start()
  s=s[:a]+('&guess_1\ng_om='+reference['initial_re']+'\ng_gam='+reference['initial_im']+'\n')+s[b:]
  (p/'case.in').write_text(s)
  (p/'preflight.txt').write_text(tensors[0]['omega_re']+' '+tensors[0]['omega_im']+'\n')
  with(p/'stdout.log').open('w')as f:
   result=subprocess.run(['mpirun','--oversubscribe','-np','4',str(root/'build/ALPS'),'case.in'],cwd=p,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=240)
  assert result.returncode==0,(p,result.returncode)
  actual=next(csv.DictReader((p/'r2_root.csv').open()))
  for key in ['omega_re','omega_im','iterations','det_re','det_im','scaled_det','ez_fraction','ey_fraction','helicity']:
   assert float(actual[key])==float(reference[key]),(kind,key,actual[key],reference[key])
  assert float(actual['scaled_det'])<=1e-11 and float(actual['ez_fraction'])<1e-5 and float(actual['ey_fraction'])<1e-8
  new=list(csv.DictReader((p/'r2_tensor.csv').open()));observed={}
  freq=(tensors[0]['omega_re'],tensors[0]['omega_im'])
  for row in new:
   if (float(row['omega_re']),float(row['omega_im']))==(float(freq[0]),float(freq[1])):
    observed[(int(row['species']),int(row['i']),int(row['j']))]=(float(row['re']),float(row['im']))
  count=0
  for row in tensors:
   if row['stage']!='1' or float(row['kpar'])!=.1:continue
   key=(int(row['species']),int(row['i']),int(row['j']))
   assert observed[key]==(float(row['re']),float(row['im'])),(kind,key)
   assert observed[(-2,0,0)]==(float(row['det_re']),float(row['det_im']))
   count+=1
  comparisons.append([kind,output,'PASS',count,actual['omega_re'],actual['omega_im'],actual['scaled_det'],actual['iterations']])
 assert (root/(kind+'_root_False')/'r2_root.csv').read_bytes()==(root/(kind+'_root_True')/'r2_root.csv').read_bytes()
 print('PASS: ordinary '+kind+' root, polarization and preflight tensors exactly match archived method-2 control',flush=True)
with(root/'root_comparison.csv').open('w')as f:
 w=csv.writer(f);w.writerow(['shape','writeOut','status','matched_susceptibility_components','omega_re','omega_im','scaled_det','iterations']);w.writerows(comparisons)

# Additional source preflight probes reuse a completed cache and the same binary.
extra=[]
for case in ['zero_degree','source_zero_cache','nonfinite_coordinate','reversed_parallel']:
 p=root/case;template=root/('cached_hit' if case=='source_zero_cache' else 'method2')
 shutil.copytree(template/'distribution',p/'distribution');(p/'solution').mkdir()
 text=(template/'case.in').read_text()
 if case=='zero_degree':text=text.replace('order=2','order=0',1)
 (p/'case.in').write_text(text)
 source=p/'distribution/susceptibility.1.array'
 if case!='zero_degree':
  lines=source.read_text().splitlines();parts=lines[0].split()
  if case=='source_zero_cache':parts[2]='0.0'
  if case=='nonfinite_coordinate':parts[0]='NaN'
  lines[0]=' '.join(parts)
  if case=='reversed_parallel':lines[0],lines[1]=lines[1],lines[0]
  source.write_text('\n'.join(lines)+'\n')
 with(p/'stdout.log').open('w')as log:
  run=subprocess.run(['mpirun','--oversubscribe','-np','4',str(root/'build/ALPS_plain'),'case.in'],cwd=p,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=60)
 log=(p/'stdout.log').read_text()
 assert run.returncode==18 and 'SIGFPE' not in log and 'ERROR: relativistic input species 1:' in log,(case,run.returncode)
 if case=='source_zero_cache':assert 'REL_SOURCE_SUPPORT' not in log
 extra.append([case,run.returncode,'PASS'])
with(root/'extra_preflight.csv').open('w')as f:
 w=csv.writer(f);w.writerow(['case','exit','status']);w.writerows(extra)
print('PASS: zero degree, invalid source with cache present, nonfinite and reversed grid checks',flush=True)

# Test-only initial exports: compare exact packed metadata/maps/value and dP/dx
# coefficient bytes on actual production ranks. Export adds no MPI operations.
signatures={}
for p in root.iterdir():
 if not p.is_dir():continue
 models=sorted(p.glob('models_rank*.bin'))
 if not models:continue
 signatures[p.name]={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in models}
 assert len(models)==4 and len(set(signatures[p.name].values()))==1,(p,signatures[p.name])
(root/'model_signatures.json').write_text(json.dumps(signatures,indent=2)+'\n')
print('PASS: production rank metadata, maps and coefficient exports byte-identical',flush=True)

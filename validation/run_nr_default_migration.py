#!/usr/bin/env python3
"""Compare production initialization with Task 3F.3 native results, no scientific sweep."""
import pathlib,sys,tempfile,subprocess,os,re,csv,json
repo=pathlib.Path(__file__).resolve().parents[1]
r=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(tempfile.mkdtemp(prefix='alps-nr-migration-'))
r=r.resolve();assert not r.is_relative_to(repo)
print('Migration artifacts:',r,flush=True)
subprocess.run([sys.executable,str(repo/'validation/build_spani_scientific.py'),str(r)],check=True)
s=(repo/'inputs/SPANi-250320-180012-v06_kpar_scan_k05c.in').read_text()
s=re.sub(r'(?im)^\s*(writeOut|fit_check|use_map)\s*=.*',lambda m:m[1]+'=F',s)
s=re.sub(r'(?im)^\s*n_scan\s*=.*','n_scan=0',s)
s=re.sub(r'(?im)^\s*(ngamma|npparbar|kperp_norm)\s*=.*','',s)
s=re.sub(r'(?im)^[ \t]*log_max[ \t]*=.*\n','',s)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
old=json.loads((repo/'audits/spani_scientific/provenance.json').read_text());reference=pathlib.Path(old['scratch'])
checks=[]
for selector in [0,1]:
 for phase in ['fixed','root']:
  p=r/f'{phase}_{selector}';(p/'distribution').mkdir(parents=True);(p/'solution').mkdir()
  (p/'case.in').write_text(s);(p/'distribution/SPANi-250320-180012-v06.1.array').symlink_to(repo/'distribution/SPANi-250320-180012-v06.1.array')
  with(p/'stdout.log').open('w')as out,(p/'stderr.log').open('w')as err:
   result=subprocess.run(['mpirun','--oversubscribe','-np','4',str(r/'build/spani'),'case.in','4','30',str(selector),phase,'1','.075'],cwd=p,env=env,stdout=out,stderr=err,timeout=240)
  if result.returncode:raise RuntimeError((p,result.returncode))
  before=reference/f'{phase}_3_N30_d{selector}_w1_k0.075_np4'
  for f in ['tensors.csv','results.csv']:
   if (before/f).exists():
    match=(p/f).read_bytes()==(before/f).read_bytes()
   else:
    with(repo/'audits/spani_scientific'/f).open()as stream:
     reader=csv.DictReader(stream);fields=reader.fieldnames[9:]
     expected=[[row[k] for k in fields] for row in reader if row['case']==before.name]
    match=list(csv.reader((p/f).open()))==expected
   checks.append([phase,selector,f,match]);assert match,(p,f)
  assert 'Determining GLLS' not in(p/'stdout.log').read_text()
with(r/'comparisons.csv').open('w')as f:w=csv.writer(f);w.writerow(['phase','selector','file','byte_identical']);w.writerows(checks)
print('PASS: production initialization reproduces native SPAN-I tensors, roots and polarization byte-identically; no GLLS fit',flush=True)
# Optional existing-fixture roots from run_nr_svd_physics.sh (Maxwellian, kappa).
# They contain unmodified generated arrays and both derivative input variants.
for fixture_root in map(pathlib.Path,sys.argv[2:]):
 for selector,kind in [(0,'nr_fd'),(1,'nr_analytic')]:
  cases=[]
  for backend in [3,4]:
   p=r/(fixture_root.name+'_'+kind+'_'+str(backend));p.mkdir();shutil.copytree(fixture_root/kind/'distribution',p/'distribution');(p/'solution').mkdir()
   (p/'case.in').write_bytes((fixture_root/kind/'case.in').read_bytes())
   with(p/'stdout.log').open('w')as out,(p/'stderr.log').open('w')as err:
    result=subprocess.run(['mpirun','--oversubscribe','-np','4',str(r/'build/spani'),'case.in',str(backend),'30',str(selector),'fixed','4','.03'],cwd=p,env=env,stdout=out,stderr=err,timeout=240)
   assert result.returncode==0,(p,result.returncode)
   assert 'Determining GLLS' not in(p/'stdout.log').read_text();cases.append(p)
  for filename in ['tensors.csv','results.csv']:
   match=(cases[0]/filename).read_bytes()==(cases[1]/filename).read_bytes();assert match
   with(r/'comparisons.csv').open('a')as f:csv.writer(f).writerow([fixture_root.name+'_'+kind,selector,filename,match])
 print('PASS: positive fixture native/production equality',fixture_root,flush=True)

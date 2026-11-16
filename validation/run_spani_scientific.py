#!/usr/bin/env python3
"""Reproducible bounded Task 3F.3 experiment. Python stdlib; Fortran does numerical work."""
import pathlib,sys,subprocess,tempfile,os,csv,json,shutil,re,math,hashlib
repo=pathlib.Path(__file__).resolve().parents[1]
root=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else pathlib.Path(tempfile.mkdtemp(prefix='alps-spani-'))
assert not root.is_relative_to(repo)
root.mkdir(parents=True,exist_ok=True);print('SPAN-I artifacts:',root,flush=True)
if not (root/'build/spani').exists():subprocess.run([sys.executable,str(repo/'validation/build_spani_scientific.py'),str(root)],check=True)
input_text=(repo/'inputs/SPANi-250320-180012-v06_kpar_scan_k05c.in').read_text()
input_text=re.sub(r'(?im)^\s*(writeOut|fit_check|use_map)\s*=.*',lambda m:m[1]+'=F',input_text)
input_text=re.sub(r'(?im)^\s*n_scan\s*=.*','n_scan=0',input_text)
input_text=re.sub(r'(?im)^\s*(ngamma|npparbar|kperp_norm)\s*=.*','',input_text)+'\n&relativistic_params\nngamma=100\nnpparbar=200\n/\n'
input_text=re.sub(r'(?im)^[ \t]*log_max[ \t]*=.*\n','',input_text)
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
manifest=[];backends=['legacy','shared','svd_index_control','native']
def run(b,n,sel,phase='fixed',iw=1,k=.075,np=4,union=None):
 name=f'{phase}_{b}_N{n}_d{sel}_w{iw}_k{k:g}_np{np}';p=root/name
 if p.exists() and (p/'exit.json').exists():
  rec=json.loads((p/'exit.json').read_text())
  if 'executable_sha256' not in rec:
   previous=root/'build/spani_initial'
   rec['executable_sha256']=hashlib.sha256((previous if previous.exists() else root/'build/spani').read_bytes()).hexdigest()
  manifest.append(rec);return p,rec
 p.mkdir(exist_ok=True);(p/'distribution').mkdir(exist_ok=True);(p/'solution').mkdir(exist_ok=True)
 dest=p/'distribution/SPANi-250320-180012-v06.1.array'
 if not dest.exists():dest.symlink_to(repo/'distribution/SPANi-250320-180012-v06.1.array')
 (p/'case.in').write_text(re.sub(r'(?im)^\s*order\s*=.*','order='+str(n),input_text))
 if union:shutil.copyfile(union,p/'query_union.csv')
 cmd=['mpirun','--oversubscribe','-np',str(np),str(root/'build/spani'),'case.in',str(b),str(n),str(sel),phase,str(iw),str(k)]
 with (p/'stdout.log').open('w') as out,(p/'stderr.log').open('w') as err:
  try:status=subprocess.run(cmd,cwd=p,env=env,stdout=out,stderr=err,timeout=240).returncode
  except subprocess.TimeoutExpired:status=124
 rec=dict(name=name,backend=backends[b],degree=n,selector=sel,phase=phase,frequency=iw,kpar=k,ranks=np,exit=status,executable_sha256=hashlib.sha256((root/'build/spani').read_bytes()).hexdigest())
 (p/'exit.json').write_text(json.dumps(rec));manifest.append(rec)
 print(name,'exit',status,flush=True);(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
 return p,rec
def result(p,label):
 try:return next(row for row in csv.reader((p/'results.csv').open()) if row[0]==label)
 except (FileNotFoundError,StopIteration):return None
# Four degrees, four fitting/evaluation pathways; native fitting diagnostics never alter the source VDF.
for n in [10,20,30,40]:
 for b in range(4):run(b,n,0,'fit',np=1)
for n in [10,20,30,40]:
 for b in range(4):
  for sel in [0,1]:
   for iw in [1,2,3]:run(b,n,sel,iw=iw)
for b in range(4):
 for sel in [0,1]:run(b,30,sel,k=.08)
# Common actual-query corpus from completed fixed-frequency cases only.
queries=set()
for m in manifest:
 if m['phase']!='fixed' or m['exit']!=0:continue
 for f in (root/m['name']).glob('queries_rank*.csv'):
  for row in csv.DictReader(f.open()):queries.add((int(row['species']),int(row['row']),float(row['u_re']),float(row['u_im'])))
union=root/'query_union.csv'
with union.open('w') as f:
 for s,r,x,y in sorted(queries):f.write(f'{s},{r},{x:.17g},{y:.17g}\n')
for n in [10,20,30,40]:
 for b in range(4):run(b,n,0,'replay',np=1,union=union)
# One documented seed per eligible configuration. No seed tuning.
roots={}
order=[(0,30,0)]+[(b,n,s) for n in [10,20,30,40] for b in range(4) for s in [0,1] if (b,n,s)!=(0,30,0)]
for b,n,s in order:
 pre=[m for m in manifest if m['phase']=='fixed' and m['kpar']==.075 and m['backend']==backends[b] and m['degree']==n and m['selector']==s]
 if len(pre)!=3 or any(m['exit']!=0 for m in pre):continue
 p,m=run(b,n,s,'root');r=result(p,'root')
 if m['exit']==0 and r:roots[b,n,s]=r
# Qualification: algebraic residual + vector character; uncertainty is reported, not hidden.
def qualified(r):return r is not None and int(r[5])<49 and float(r[8])<=1e-7
def overlap(a,b):
 ea=[complex(float(a[9+i]),float(a[12+i])) for i in range(3)];eb=[complex(float(b[9+i]),float(b[12+i])) for i in range(3)]
 return abs(sum(x.conjugate()*y for x,y in zip(ea,eb)))/math.sqrt(sum(abs(x)**2 for x in ea)*sum(abs(y)**2 for y in eb))
for s in [0,1]:
 ref=roots.get((0,30,s));nat=roots.get((3,30,s))
 if qualified(ref) and qualified(nat) and overlap(ref,nat)>=.95:
  for b,n in [(0,30),(1,30),(2,30),(3,30),(3,10),(3,20),(3,40)]:
   if qualified(roots.get((b,n,s))):run(b,n,s,'scan')
# Representative actual assembly comparisons, with both derivative selectors.
for b in [0,3]:
 for s in [0,1]:
  for np in [1,2]:
   run(b,30,s,np=np)
   if qualified(roots.get((b,30,s))):run(b,30,s,'root',np=np)
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('COMPLETE:',len(manifest),'case records;',len(queries),'unique actual queries',flush=True)

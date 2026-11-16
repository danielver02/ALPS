#!/usr/bin/env python3
"""Ordinary production-input error checks; private copies of the SPAN-I table."""
import pathlib,sys,shutil,subprocess,os,tempfile,re
exe=pathlib.Path(sys.argv[1]).resolve();source=pathlib.Path(sys.argv[2]).resolve();root=pathlib.Path(tempfile.mkdtemp(prefix='alps-nr-default-errors-'))
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for kind,code in [('negative',16),('unavailable',15)]:
 p=root/kind;(p/'distribution').mkdir(parents=True);(p/'solution').mkdir()
 (p/'case.in').write_text(re.sub(r'(?im)^\s*writeOut\s*=.*','writeOut=F',(source/'case.in').read_text()))
 src=next((source/'distribution').glob('*.1.array'));lines=src.read_text().splitlines()
 for j in ([402] if kind=='negative' else range(402,804)):
  a=lines[j].split();a[2]='-0.1' if kind=='negative' else '0';lines[j]=' '.join(a)
 (p/'distribution'/src.name).write_text('\n'.join(lines)+'\n')
 with(p/'stdout.log').open('w')as out,(p/'stderr.log').open('w')as err:
  r=subprocess.run(['mpirun','--oversubscribe','-np','4',str(exe),'case.in'],cwd=p,env=env,stdout=out,stderr=err,timeout=120)
 text=(p/'stdout.log').read_text();assert r.returncode!=0 and re.search(r'iperp=1 status='+str(code)+r'\b',text),(p,r.returncode)
 assert 'Determining GLLS'not in text
 assert 'SIGFPE'not in(p/'stderr.log').read_text()
 print('PASS ordinary production:',kind,'status',code,flush=True)
print('Error artifacts:',root)

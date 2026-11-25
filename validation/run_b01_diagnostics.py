#!/usr/bin/env python3
"""Scratch ordinary-production writeOut check. Arguments: output, SPAN-I case, NR case."""
import pathlib,sys,subprocess,shutil,re,os,csv,math
repo=pathlib.Path(__file__).resolve().parents[1]
out=pathlib.Path(sys.argv[1]).resolve();assert not out.is_relative_to(repo)
out.mkdir();(out/'source').mkdir();(out/'build').mkdir()
for p in (repo/'src').glob('*.f90'):shutil.copy2(p,out/'source'/p.name)
shutil.copy2(repo/'distribution/distribution_analyt.f90',out/'source')
names=['ALPS_rel_cache','ALPS_var','distribution_analyt','ALPS_io','ALPS_analyt','ALPS_com','ALPS_NHDS','ALPS_fns_rel','ALPS_fns','ALPS_nyquist','ALPS_nyquist_growth','ALPS_nyquist_driver']
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
rows=[]
for opt in ['-O0','-O4']:
 build=out/'build'
 flags=[opt,'-g','-fcheck=bounds','-ffree-line-length-none','-ffpe-trap=overflow']
 for name in names:subprocess.run(['mpifort',*flags,'-c',str(out/'source'/(name+'.f90'))],cwd=build,check=True)
 subprocess.run(['mpifort',*flags,'-o','ALPS',*[n+'.o' for n in names],str(out/'source/ALPS.f90'),'-llapack','-lblas'],cwd=build,check=True)
 for kind,fixture in [('native',sys.argv[2]),('functional',sys.argv[3]),('analytic',sys.argv[3])]:
  fixture=pathlib.Path(fixture).resolve();outputs=[]
  for enabled in [False,True]:
   case=out/(opt[1:]+'_'+kind+'_'+str(enabled));case.mkdir();(case/'solution').mkdir()
   shutil.copytree(fixture/'distribution',case/'distribution')
   for p in (case/'distribution').glob('*.out'):p.unlink()
   text=(fixture/'case.in').read_text()
   text=re.sub(r'(?im)^[ \t]*log_max[ \t]*=.*\n','',text)
   text=re.sub(r'(?im)^\s*writeOut\s*=.*','writeOut='+('T' if enabled else 'F'),text)
   if kind!='native':text=re.sub(r'(?im)^\s*AC_method\s*=.*','AC_method='+('1' if kind=='functional' else '0'),text)
   if kind!='native':
    text=re.sub(r'(?im)^\s*use_map\s*=.*','use_map=T',text)
    for key,value in [('omi','1.55d0'),('omf','1.56d0'),('gami','.01d0'),('gamf','.02d0'),('ni','2'),('nr','2')]:
     text=re.sub(r'(?im)^\s*'+key+r'\s*=.*',key+'='+value,text)
    # A separate small analytic Maxwellian diagnostic fixture; not a resampling
    # or alteration of any archived numerical reference or measured VDF.
    text=re.sub(r'(?im)^\s*nperp\s*=.*','nperp=8',text)
    text=re.sub(r'(?im)^\s*npar\s*=.*','npar=32',text)
    for species,width,amplitude in [(1,1.,math.pi**(-1.5)),(2,1./math.sqrt(1836.),(1836./math.pi)**1.5)]:
     with(case/'distribution'/('test_bimax.'+str(species)+'.array')).open('w')as stream:
      for row in range(9):
       vp=3.*width*row/8.
       for node in range(33):
        vu=3.*width*(-1.+node/16.)
        stream.write(f'{vp:.17e} {vu:.17e} {amplitude*math.exp(-(vp*vp+vu*vu)/(width*width)):.17e}\n')

   (case/'case.in').write_text(text)
   with (case/'stdout.log').open('w')as stdout,(case/'stderr.log').open('w')as stderr:
    run=subprocess.run(['mpirun','--oversubscribe','-np','4',str(build/'ALPS'),'case.in'],cwd=case,env=env,stdout=stdout,stderr=stderr,timeout=240)
   assert run.returncode==0,(case,run.returncode)
   log=(case/'stdout.log').read_text()
   assert 'Determining GLLS'not in log
   if enabled:
    if kind=='native':
     assert not list((case/'distribution').glob('*.poly_fit_result.*'))
     assert 'Native SVD rows excluded' in log and 'Bi-Maxwellian/cold-plasma Species' in log
    else:assert list((case/'distribution').glob('*.fit_result.*'))
   outputs.append({p.name:p.read_bytes() for p in (case/'solution').iterdir()if p.is_file()})
  assert outputs[0]==outputs[1],(opt,kind,'scientific output changed with writeOut')
  rows.append([opt,kind,'PASS',len(outputs[0])])
with(out/'results.csv').open('w')as f:w=csv.writer(f);w.writerow(['optimization','case','status','identical_solution_files']);w.writerows(rows)
print(rows,flush=True)

#!/usr/bin/env python3
"""Exercise actual init_param in isolated MPI jobs, including unused-method blocks."""
from pathlib import Path
import csv, os, shlex, subprocess, tempfile
from lm_input_fixture import relocate_lm_controls
repo=Path(__file__).resolve().parents[1]
root=Path(tempfile.mkdtemp(prefix='alps-fitting-input-'))
print(f'Fitting input artifacts: {root}',flush=True)
fc=shlex.split(os.environ.get('MPIFC','mpifort'))
flags=shlex.split(os.environ.get('FFLAGS','-O0 -g -fcheck=all -fbacktrace'))
subprocess.run(fc+flags+['-o','input',str(repo/'src/ALPS_var.f90'),str(repo/'src/ALPS_io.f90'),str(repo/'validation/test_fitting_input.f90')],cwd=root,check=True)
mpi=shlex.split(os.environ.get('MPIEXEC','mpirun'))+shlex.split(os.environ.get('MPIEXEC_FLAGS','--oversubscribe'))
system='&system nspec=2,nroots=0,use_map=F,n_scan=0,nperp=4,npar=8,writeOut=F /\n'
explicit='&fitting_params maxsteps_fit=500,lambda_initial_fit=1.d0,lambdafac_fit=1.d1,epsilon_fit=1.d-8 /\n'
nondefault='&fitting_params maxsteps_fit=17,lambda_initial_fit=0.25d0,lambdafac_fit=3.d0,epsilon_fit=2.d-9 /\n'
def species(methods,rel=False):
 s=''
 for j,m in enumerate(methods,1):
  s+=f'&spec_{j}\n nn=1,qq=1,mm=1,ff=1,AC_method={m},relat={"T" if rel else "F"},log_fit=T,use_bM=F /\n'
  if m==1:s+=f'&ffit_{j}_1\n fit_type_in={4 if rel else 1},fit_1=1,fit_2=1,fit_3=0,fit_4=0,fit_5=0,perpcorr=2 /\n'
  if m==2:s+=f'&poly_spec_{j}\n kind=1,order=2,derivative_type=0 /\n'
 return s
cases=[]
def add(name,text,expected=None,error=None,methods=(1,1),grid=(100,200,0)):
 cases.append((name,text,expected,error,methods,grid))
default=(500,1.,10.,1.e-8);custom=(17,.25,3.,2.e-9)
for methods in [(0,0),(1,1),(2,2),(1,2)]:
 label=''.join(map(str,methods));tail=species(methods)
 add('omitted_'+label,system+tail,default,methods=methods)
 add('explicit_'+label,system+explicit+tail,default,methods=methods)
 add('nondefault_'+label,system+nondefault+tail,custom,methods=methods)
add('partial',system+'&fitting_params maxsteps_fit=17 /\n'+species((1,1)),(17,1.,10.,1.e-8))
add('rel_mixed',system+species((1,2),True)+nondefault+'&relativistic_params ngamma=8,npparbar=12 /\n',custom,methods=(1,2),grid=(8,12,0))
add('before_system',explicit+system+species((1,1)),default)
add('after_optional',system+species((1,2),True)+'&relativistic_params ngamma=8,npparbar=12 /\n'+nondefault,custom,methods=(1,2),grid=(8,12,0))
add('comments_strings',system.replace('writeOut=F','arrayName="text &fitting_params",writeOut=F')+'! &fitting_params /\n'+species((1,1))+"! quote '&fitting_params'\n",default)
add('case_inline',system+'   &FiTtInG_PaRaMs maxsteps_fit=17 /\n'+species((1,1)),(17,1.,10.,1.e-8))
add('long_comment',system+'! '+('x'*3000)+' &fitting_params /\n'+explicit+species((1,1)),default)
add('quoted_escape',system.replace('writeOut=F',"arrayName='it''s &fitting_params',writeOut=F")+explicit+species((1,1)),default)
# Numeric range policy is deliberately unchanged. Parse only; do not run LM with invalid controls.
add('empty',system+'&fitting_params /\n'+species((1,1)),default)
scan='&scan_input_1\n scan_type=4,swi=.1,swf=.2,swlog=F,ns=2,nres=1,heating=F,eigen=F /\n'
for position in ['before','after']:
 tail=species((1,2),True)+'&relativistic_params ngamma=8,npparbar=12 /\n'
 tail=(nondefault+tail+scan) if position=='before' else (tail+scan+nondefault)
 add('scan_'+position,system.replace('n_scan=0','n_scan=1')+tail,custom,methods=(1,2),grid=(8,12,1))
add('range_policy',system+'&fitting_params maxsteps_fit=-1,lambda_initial_fit=-2,lambdafac_fit=0,epsilon_fit=-1 /\n'+species((1,1)),(-1,-2.,0.,-1.))
legacy_lines='\n'.join(f'{key}={value}' for key,value in [('maxsteps_fit',17),('lambda_initial_fit','.25d0'),('lambdafac_fit','3.d0'),('epsilon_fit','2.d-9')])
legacy=system.replace('&system ','&system\n ').replace(' /','\n'+legacy_lines+'\n/')+species((1,1))
legacy=legacy.replace('\n/','\nfit_check=F\n/',1)
converted=relocate_lm_controls(legacy)
assert 'fit_check' not in converted
assert relocate_lm_controls(converted)==converted
add('converted_archived_fixture',converted,custom)
try:
 relocate_lm_controls(legacy+explicit)
except ValueError:
 pass
else:
 raise AssertionError('ambiguous archived fixture conversion accepted')
for key,value in [('maxsteps_fit','500'),('lambda_initial_fit','1.d0'),('lambdafac_fit','10.d0'),('epsilon_fit','1.d-8'),('fit_check','T')]:
 add('legacy_'+key,system.replace(' /',f',{key}={value} /')+explicit+species((1,1)),error='to &fitting_params')
for name,block,error in [
 ('duplicate',explicit+nondefault,'Duplicate &fitting_params'),
 ('unknown','&fitting_params not_a_control=1 /\n','Cannot read &fitting_params'),
 ('bad_integer','&fitting_params maxsteps_fit="bad" /\n','Cannot read &fitting_params'),
 ('bad_real','&fitting_params epsilon_fit="bad" /\n','Cannot read &fitting_params'),
 ('removed_setting','&fitting_params fit_check=T /\n','Cannot read &fitting_params'),
 ('truncated','&fitting_params maxsteps_fit=17\n','Cannot read &fitting_params'),
 ('indexed','&fitting_params_1 maxsteps_fit=17 /\n','without a species index'),
]:add(name,system+species((1,1))+block,error=error)
add('nested',system+species((1,1)).replace(' nn=1',' &fitting_params maxsteps_fit=17 /\n nn=1',1),error='between complete namelist blocks')
for methods in [(0,0),(2,2)]:
 add('invalid_unused_'+str(methods[0]),system+species(methods)+'&fitting_params bogus=1 /\n',error='Cannot read &fitting_params',methods=methods)
# An omitted per-species selector uses the analytic default independently of
# explicit settings in preceding species blocks.
poly_pair=species((2,2))
add('poly_default',system+poly_pair.replace(',derivative_type=0',''),default,methods=(2,2))
add('poly_explicit_analytic',system+poly_pair.replace('derivative_type=0','derivative_type=1'),default,methods=(2,2))
add('poly_mixed_default',system+poly_pair.replace(',derivative_type=0','',1),default,methods=(2,2))
rows=[]
for name,text,expected,error,methods,grid in cases:
 p=root/name;p.mkdir();(p/'case.in').write_text(text)
 with (p/'stdout.log').open('w') as out,(p/'stderr.log').open('w') as err:
  result=subprocess.run(mpi+['-np','4',str(root/'input'),'case.in'],cwd=p,stdout=out,stderr=err,timeout=60)
 stdout=(p/'stdout.log').read_text();stderr=(p/'stderr.log').read_text()
 if error:
  assert result.returncode!=0,(name,'accepted invalid input')
  assert error in stdout,(name,stdout,stderr)
  assert 'Fortran runtime error' not in stdout+stderr,(name,stderr)
 else:
  assert result.returncode==0,(name,stdout,stderr)
  assert 'PASS: root-only' in stdout,(name,stdout)
  lines=(p/'parsed.dat').read_text().splitlines();data=lines[0].split()
  values=(int(data[0]),*map(lambda s:float(s.replace('D','E')),data[1:4]))
  assert values==expected,(name,values,expected)
  assert tuple(map(int,lines[1].split()))==methods,name
  assert tuple(map(int,lines[2].split()))==tuple(2 if m==2 else 0 for m in methods),name
  assert tuple(map(int,lines[3].split()))==grid,name
  selectors=(0,0)  # Non-polynomial species retain their finite-difference path.
  if name in ['poly_default','poly_explicit_analytic']:selectors=(1,1)
  if name=='poly_mixed_default':selectors=(1,0)
  assert tuple(map(int,lines[4].split()))==selectors,(name,lines[4],selectors)
 rows.append([name,result.returncode,'PASS']);print('PASS:',name,flush=True)
with(root/'results.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['case','exit','result']);w.writerows(rows)
print(f'PASS: {len(rows)} real-parser MPI cases; artifacts: {root}',flush=True)

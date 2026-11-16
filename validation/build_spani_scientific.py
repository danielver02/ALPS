#!/usr/bin/env python3
"""Build copied ALPS sources; instrument only the copies. Standard library only."""
import pathlib,subprocess,shutil,re,sys
repo=pathlib.Path(__file__).resolve().parents[1];out=pathlib.Path(sys.argv[1]).resolve()
assert not out.is_relative_to(repo)
(out/'source').mkdir(parents=True);(out/'build').mkdir()
for p in (repo/'src').glob('*.f90'):shutil.copy2(p,out/'source'/p.name)
for p in [repo/'distribution/distribution_analyt.f90',repo/'validation/legacy_glls_reference.f90',repo/'validation/spani_diagnostics.f90',repo/'validation/test_spani_scientific.f90']:shutil.copy2(p,out/'source'/p.name)
subprocess.run([sys.executable,str(repo/'validation/instrument_nr_history.py'),str(out/'source')],check=True)
p=out/'source/ALPS_analyt.f90';s=p.read_text()
a=' eval_fit=historical_value(is,iperp,ppar_valC,present(dfpar),historical_df)';assert s.count(a)==1;s=s.replace(a,a+'\n call span_query(is,iperp,ppar_valC,present(dfpar),0)')
a='   call eval_nr_cheb(is,iperp,ppar_valC,eval_fit,cheb_status,dfpar,value_status,derivative_status)'
assert s.count(a)==1;s=s.replace(a,a+'\n   call span_query(is,iperp,ppar_valC,present(dfpar),cheb_status)')
a='end select\nend function eval_fit';assert s.count(a)==1;s=s.replace(a,'end select\ncall span_query(is,iperp,ppar_valC,present(dfpar),cheb_status)\nend function eval_fit');p.write_text(s)
p=out/'source/ALPS_fns.f90';s=p.read_text();start=s.index(' double complex function disp(om)');end=s.index('end function disp',start)+len('end function disp');serial=s[start:end]
serial=serial.replace('function disp','function span_serial_disp').replace('disp =','span_serial_disp =').replace('mpi_bcast(disp,','mpi_bcast(span_serial_disp,')
serial=serial.replace('only : usebM','only : usebM,nmax')
a='    else\n        ! Integrate:';assert serial.count(a)==1
serial=serial.replace(a,'    endif\n    do sproc=1,nspec\n      nlim=[0,nmax(sproc)]\n        ! Integrate:')
a='    endif\n\n    ! Return the schi';assert serial.count(a)==1;serial=serial.replace(a,'    enddo\n\n    ! Return the schi')
serial,n=re.subn(r'    call MPI_REDUCE \(schi, chi, size\(chi\),&\s*MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror\)','    chi=schi',serial);assert n==1
serial,n=re.subn(r'    call MPI_REDUCE \(schi_low, chi_low, size\(chi_low\),&\s*MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierror\)','    chi_low=schi_low',serial);assert n==1
serial=serial.replace('    return\n','    sproc=0\n    return\n')
s=s.replace('end module alps_fns',serial+'\nend module alps_fns') if 'end module alps_fns' in s else s.replace('end module ALPS_fns',serial+'\nend module ALPS_fns')
assert s.count('function span_serial_disp')==2
s=s.replace('        f=func(rtsec)','        f=func(rtsec)\n        call span_step(j)')
s=s.replace('           tmp = disp(omega)','           tmp = disp(omega)\n           call span_scan_result(omega,tmp)',1)
p.write_text(s)
fc='mpifort';flags=['-O2','-g','-fcheck=bounds','-fbacktrace','-ffree-line-length-none','-ffpe-trap=invalid,zero,overflow']
names=['ALPS_rel_cache','ALPS_var','distribution_analyt','ALPS_io','ALPS_analyt','ALPS_com','legacy_glls_reference','ALPS_NHDS','ALPS_fns_rel','ALPS_fns','spani_diagnostics']
# Diagnostics depends on ALPS_fns, whereas scratch hooks use external procedures.
for n in names:subprocess.run([fc,*flags,'-c','../source/'+n+'.f90'],cwd=out/'build',check=True)
subprocess.run([fc,*flags,'-o','spani',*[n+'.o' for n in names],'../source/test_spani_scientific.f90','-llapack','-lblas'],cwd=out/'build',check=True)

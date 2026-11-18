#!/usr/bin/env python3
"""Install historical dispatch only in a copied validation source tree."""
import pathlib,sys
root=pathlib.Path(sys.argv[1]).resolve();repo=pathlib.Path(__file__).resolve().parents[1]
assert not root.is_relative_to(repo),'instrumentation must remain outside repository'
p=next(root.rglob('ALPS_analyt.f90'));s=p.read_text();a=s.index('double complex function eval_fit(');b=s.index('end function eval_fit',a)
part=s[a:b];part=part.replace('implicit none','implicit none\n logical,external::historical_active\n double complex,external::historical_value\n double complex::historical_df',1)
anchor='if(ACmethod(is)==2.and.relativistic(is)) then';assert part.count(anchor)==1
part=part.replace(anchor,'''if(historical_active(is))then
 eval_fit=historical_value(is,iperp,ppar_valC,present(dfpar),historical_df)
 if(present(dfpar))dfpar=historical_df
 return
endif
'''+anchor)
s=s[:a]+part+s[b:];p.write_text(s)
p=next(root.rglob('ALPS_com.f90'));s=p.read_text();anchor='  end subroutine pass_distribution';assert s.count(anchor)==1;s=s.replace(anchor,'    call historical_distribution ! Validation-owned transfer only.\n'+anchor);p.write_text(s)

# Only validation executables retain the historical method-2 contour gate.
# Instrument the actual routine, not the generated before/candidate or physics.
p=next(root.rglob('ALPS_fns.f90'));s=p.read_text()
a=s.index('double complex function landau_integrate(');b=s.index('end function landau_integrate',a)
part=s[a:b];part=part.replace('implicit none','implicit none\n logical,external::historical_active',1)
anchor='legacy_zero_return=(ACmethod(sproc)==0.or.ACmethod(sproc)==1)'
assert part.count(anchor)==1
part=part.replace(anchor,anchor+'.or.historical_active(sproc)')
s=s[:a]+part+s[b:];p.write_text(s)

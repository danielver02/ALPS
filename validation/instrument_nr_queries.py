#!/usr/bin/env python3
"""Trace/inject queries only in copied zero-contour validation sources."""
import pathlib,re,sys
root=pathlib.Path(sys.argv[1]).resolve();repo=pathlib.Path(__file__).resolve().parents[1]
assert not root.is_relative_to(repo)
p=root/'ALPS_analyt.f90';s=p.read_text();a=s.index('double complex function eval_fit(');b=s.index('end function eval_fit',a)
part=s[a:b].replace('implicit none','implicit none\n logical::b5_injected',1)
anchor='if(present(status)) status=CHEB_SUCCESS'
assert part.count(anchor)==1
part=part.replace(anchor,"""call b5_query(is,iperp,ppar_valC,present(dfpar),b5_injected)
if(b5_injected)call alps_continuation_error(is,iperp,CHEB_LOG_OVERFLOW,ppar_valC,nonrelativistic=.true.)
"""+anchor)
s=s[:a]+part+s[b:];p.write_text(s)
p=root/'ALPS_fns.f90';s=p.read_text()
for name in ['landau_integrate','landau_integrate_before','landau_integrate_candidate']:
 a=s.index('double complex function '+name+'(');b=s.index('end function '+name,a)
 part=s[a:b]
 part,count=re.subn(r'(^[ \t]*p_res = .*?$)',r'\1\n  call b5_row(iperp)',part,flags=re.M)
 assert count==3,(name,count)
 s=s[:a]+part+s[b:]
p.write_text(s)

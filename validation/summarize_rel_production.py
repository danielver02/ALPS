#!/usr/bin/env python3
"""Independent finite-grid quadrature/shape diagnostics on read-only R2 dumps."""
from pathlib import Path
import csv,math,sys
root=Path(sys.argv[1]).resolve();out=root/'analysis';out.mkdir(exist_ok=True)
metrics=[]
for case,b in [('method2',.15),('iso_root_False',0.),('aniso_root_False',.15)]:
 rows=list(csv.DictReader((root/case/'r1_nodes.csv').open()))
 for bounds in csv.DictReader((root/case/'r1_bounds.csv').open()):
  species=bounds['species'];data=[r for r in rows if r['species']==species]
  ng=int(bounds['ngamma']);nu=int(bounds['npparbar']);mass=float(bounds['mass']);va=float(bounds['va'])
  dg=(float(bounds['gamma_max'])-float(bounds['gamma_min']))/ng
  du=(float(bounds['ppar_max'])-float(bounds['ppar_min']))*va/mass/nu
  def shape(r):return math.exp(-2*float(r['gamma'])-b*float(r['u'])**2)
  def integral(value):
   total=0.
   for g in range(ng+1):
    row=sorted([r for r in data if int(r['g'])==g],key=lambda r:int(r['j']))
    # Same declared finite-grid quadrature; derive from node coordinates,
    # independently of the constructor's pre-normalization accumulator.
    for left,right in zip(row,row[1:]):
     assert int(right['j'])==int(left['j'])+1
     width=float(right['u'])-float(left['u'])
     total+=(.5 if g in (0,ng) else 1.)*dg*.5*float(left['gamma'])*width*(value(left)+value(right))
   return total*2*math.pi*(mass/va)**3
  norm=integral(shape);density=integral(lambda r:float(r['table_value']))
  peak=max(shape(r)/norm for r in data)
  table_errors=[];inside=[];outside=[];fit=[];fu=[]
  for r in data:
   assert int(r['fit_status'])==0
   exact=shape(r)/norm;err=abs(float(r['table_value'])-exact)/exact
   table_errors.append(err);(inside if r['inside_source']=='T' else outside).append(err)
   fit.append(abs(float(r['fit_value'])-float(r['table_value']))/float(r['table_value']))
   fu.append(abs(float(r['fit_derivative'])+2*b*float(r['u'])*exact)/peak)
  metrics.append([case,species,len(data),len(outside),norm,density,max(table_errors),max(inside),max(outside),max(fit),max(fu)])
with(out/'analytic_shapes.csv').open('w')as f:
 w=csv.writer(f);w.writerow(['case','species','cone_nodes','outside_source','analytic_shape_integral','completed_discrete_density','max_table_vs_shape_relative','max_source_inside_relative','max_source_exterior_relative','max_fit_vs_table_relative','max_Fu_error_over_shape_peak']);w.writerows(metrics)
print('PASS: independent finite-grid shape diagnostics',flush=True)
assert all(abs(row[5]-1)<2e-12 for row in metrics)

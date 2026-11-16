#!/usr/bin/env python3
"""Aggregate completed SPAN-I scratch outputs; no reference updates. Stdlib + gnuplot."""
import pathlib,sys,csv,json,math,gzip,subprocess,hashlib,collections
root=pathlib.Path(sys.argv[1]).resolve();out=pathlib.Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
repo=pathlib.Path(__file__).resolve().parents[1]
manifest=json.loads((root/'manifest.json').read_text());(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
raw=[]
for l in (repo/'distribution/SPANi-250320-180012-v06.1.array').open():raw.append(tuple(map(float,l.split())))
perps=sorted(set(a for a,b,c in raw));us=[b for a,b,c in raw if a==perps[0]];m=len(us);F=[raw[i*m:(i+1)*m] for i in range(len(perps))]
def weights(v):return [.5*(v[min(i+1,len(v)-1)]-v[max(i-1,0)]) for i in range(len(v))]
wp,wu=weights(perps),weights(us)
rho=[sum(f*wu[j] for j,(p,u,f) in enumerate(row))*2*math.pi*perps[i]*wp[i] for i,row in enumerate(F)]
cdf=[];s=0
for v in rho:s+=v;cdf.append(s)
rows={next(i for i,v in enumerate(cdf) if v>=q*cdf[-1]) for q in [.05,.25,.5,.75,.95]}
curv=[rho[i]*max(abs(math.log10(row[j+1][2])-2*math.log10(row[j][2])+math.log10(row[j-1][2])) for j in range(1,m-1)) for i,row in enumerate(F)]
rows.add(max(range(len(curv)),key=curv.__getitem__))
counts=collections.Counter('nonfinite' if not all(map(math.isfinite,x)) else 'positive' if x[2]>0 else 'zero' if x[2]==0 else 'negative' for x in raw)
qual=dict(points=len(raw),perpendicular_points=len(perps),parallel_points=m,perpendicular_range=[perps[0],perps[-1]],parallel_range=[us[0],us[-1]],counts=counts,minimum=min(x[2] for x in raw),maximum=max(x[2] for x in raw),parallel_step_range=[min(b-a for a,b in zip(us,us[1:])),max(b-a for a,b in zip(us,us[1:]))],uniform_deviation=max(abs(u-(us[0]+(us[-1]-us[0])*j/(m-1))) for j,u in enumerate(us)),representative_rows=sorted(rows),density_trapezoidal=sum(rho))
(out/'qualification.json').write_text(json.dumps(qual,indent=2))
with (out/'vdf.dat').open('w') as f:
 for row in F:
  for p,u,v in row:f.write(f'{u:.17g} {p:.17g} {math.log10(v):.17g}\n')
  f.write('\n')
fit_summary=[];coeff=[];cuts=[];moments=[];rowdata=[]
for rec in manifest:
 if rec['phase']!='fit' or rec['exit']!=0:continue
 p=root/rec['name'];pre=[rec['backend'],rec['degree']]
 for r in csv.DictReader((p/'fit_rows.csv').open()):rowdata.append(dict(backend=pre[0],degree=pre[1],**r))
 for r in csv.DictReader((p/'coefficients.csv').open()):
  if int(r['row']) in rows:coeff.append([*pre,r['row'],r['order'],r['coefficient']])
 for r in csv.DictReader((p/'moments.csv').open()):
  d,a,b,c=map(float,[r['density'],r['first_parallel'],r['second_parallel'],r['second_perpendicular']]);moments.append([*pre,r['kind'],d,a/d,b/d-(a/d)**2,c/(2*d)])
 L1=L2=logsum=0;bad=0;maxabs=0;peak=qual['maximum']
 for r in csv.DictReader((p/'nodes.csv').open()):
  i,j=int(r['row']),int(r['node']);v,w=map(float,[r['original'],r['reconstructed']]);weight=wp[i]*wu[j]*2*math.pi*perps[i]
  if int(r['status'])!=0 or not math.isfinite(w):bad+=1;continue
  L1+=weight*abs(w-v);L2+=weight*(w-v)**2;maxabs=max(maxabs,abs(w-v))
  if w>0:logsum+=weight*v*(math.log10(w)-math.log10(v))**2
  else:bad+=1
  if i in rows:cuts.append([*pre,i,j,perps[i],us[j],v,w])
 fit_summary.append([*pre,L1/sum(rho),math.sqrt(L2),maxabs/peak,math.sqrt(logsum/sum(rho)),bad])
def write(name,header,rows):
 with (out/name).open('w') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)
write('fit_summary.csv',['backend','degree','weighted_relative_L1','absolute_L2','max_error_over_global_peak','density_weighted_log_RMS','invalid_or_zero_nodes'],fit_summary)
write('moments.csv',['backend','degree','kind','density','mean_u','variance_u','perpendicular_second_per_component'],moments)
write('coefficients.csv',['backend','degree','row','order','coefficient'],coeff)
write('cuts.csv',['backend','degree','row','node','pperp','u','original','reconstructed'],cuts)
if rowdata:write('rows.csv',list(rowdata[0]),[list(r.values()) for r in rowdata])
results=[];tensors=[];querydata=[]
for rec in manifest:
 p=root/rec['name'];prefix=[rec[k] for k in ['name','backend','degree','selector','phase','frequency','kpar','ranks','exit']]
 if (p/'results.csv').exists():
  for row in csv.reader((p/'results.csv').open()):results.append(prefix+row)
 if (p/'tensors.csv').exists():
  for row in csv.reader((p/'tensors.csv').open()):tensors.append(prefix+row)
 if rec['phase']=='replay' and(p/'replay.csv').exists():
  # Preserve the complete common corpus, compressed without any filtering.
  with gzip.open(out/(rec['backend']+'_N'+str(rec['degree'])+'_queries.csv.gz'),'wb') as f:f.write((p/'replay.csv').read_bytes())
 if (p/'stdout.log').exists():
  text=(p/'stdout.log').read_text()
  if 'ERROR:' in text or rec['exit']!=0:querydata.append([rec['name'],rec['exit'],' | '.join(l for l in text.splitlines() if 'ERROR:' in l or 'SPAN_BRANCH_STOP' in l)])
prefix=['case','backend','degree','selector','phase','probe','initial_kpar','ranks','exit']
write('results.csv',prefix+['label','point','kpar','omega_re','omega_im','iterations','det_re','det_im','scaled_residual','Ex_re','Ey_re','Ez_re','Ex_im','Ey_im','Ez_im','Bx_re','By_re','Bz_re','Bx_im','By_im','Bz_im'],results)
write('tensors.csv',prefix+['label','point','kpar','omega_re','omega_im','species','i','j','value_re','value_im'],tensors)
write('failures.csv',['case','exit','diagnostic'],querydata)
# Matched fixed-frequency component comparisons, tensor-scaled using the legacy species tensor.
lookup={}
for r in tensors:
 if r[4]!='fixed' or r[8]!=0:continue
 key=(r[2],r[3],r[5],r[6],r[7],r[14],r[15],r[16]);lookup[r[1],key]=complex(float(r[17]),float(r[18]))
comparison=[]
for (backend,key),v in lookup.items():
 if ('legacy',key)not in lookup:continue
 ref=lookup['legacy',key];d=v-ref
 family=key[:5]
 scale=max([abs(w) for (b,k),w in lookup.items() if b=='legacy' and k[:5]==family and (int(k[5])==0 if int(key[5])==0 else int(k[5])>0)]+[1e-300])
 comparison.append([backend,*key,v.real,v.imag,ref.real,ref.imag,d.real,d.imag,d.real/scale,d.imag/scale])
write('tensor_differences.csv',['backend','degree','selector','probe','kpar','ranks','species','i','j','candidate_re','candidate_im','legacy_re','legacy_im','delta_re','delta_im','tensor_scaled_re','tensor_scaled_im'],comparison)
# Gnuplot inputs: include original cuts once, and native degree curves for each selected row.
for row in sorted(rows):
 with (out/f'cut_{row}.dat').open('w') as f:
  a={n:{int(r[3]):float(r[7]) for r in cuts if r[0]=='native' and r[1]==n and r[2]==row} for n in [10,20,30,40]}
  for j,u in enumerate(us):f.write(' '.join(map(str,[u,F[row][j][2],*[a[n].get(j,float('nan')) for n in a]]))+'\n')
with (out/'spectra.dat').open('w') as f:
 for r in coeff:
  if r[0]=='native':f.write(' '.join(map(str,[r[1],r[2],r[3],abs(float(r[4]))]))+'\n')
plot="set terminal svg size 1150,750 dynamic\nset output 'vdf_cuts.svg'\nset multiplot layout 2,3\nset xlabel 'parallel momentum';set ylabel 'F';set logscale y\n"
for row in sorted(rows):plot+=f"set title 'perpendicular row {row}';plot 'cut_{row}.dat' u 1:2 w l t 'table',for [c=3:6] '' u 1:c w l t sprintf('native N%d',10*(c-2))\n"
plot+="unset multiplot\nunset logscale y\nset output 'vdf_overview.svg'\nset view map;set pm3d map;set xlabel 'parallel momentum';set ylabel 'perpendicular momentum';set cblabel 'log10 F'\nsplot 'vdf.dat' u 1:2:3 notitle\nunset pm3d\nset output 'coefficients.svg'\nset logscale y;set xlabel 'Chebyshev order';set ylabel '|coefficient|'\nplot for [n=10:40:10] 'spectra.dat' u ($1==n?$3:1/0):4 w p t sprintf('N%d (selected rows)',n)\n"
(out/'figures.gp').write_text(plot);subprocess.run(['gnuplot','figures.gp'],cwd=out,check=True)
(out/'provenance.json').write_text(json.dumps(dict(scratch=str(root),HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),array_sha256=hashlib.sha256((repo/'distribution/SPANi-250320-180012-v06.1.array').read_bytes()).hexdigest(),source_sha256={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (repo/'src').glob('*.f90')}),indent=2))
print('Packaged',len(manifest),'cases;',len(results),'result records;',len(tensors),'tensor records',flush=True)
# Common-query diagnostics: highest tested degree is a comparison, never ground truth.
replays={}
for n in [10,20,30,40]:
 for b in ['legacy','shared','svd_index_control','native']:
  p=out/(b+'_N'+str(n)+'_queries.csv.gz')
  if p.exists():
   with gzip.open(p,'rt') as f:replays[b,n]={tuple(r[k] for k in ['species','row','u_re','u_im']):r for r in csv.DictReader(f)}
summary=[];behavior=[];terms=[]
def complex_field(r,tag):return complex(float(r[tag+'_re']),float(r[tag+'_im']))
if ('native',40) in replays:
 refq=replays['native',40]
 selected={}
 for k,r in refq.items():
  row=int(r['row']);score=(-abs(float(r['u_im'])),abs(float(r['u_re'])+.067722/.075))
  if row not in selected or score<selected[row][0]:selected[row]=(score,k)
 for (b,n),items in replays.items():
  statuses=collections.Counter(r['status'] for r in items.values());outside=sum(r['outside_real']=='T' for r in items.values());maxv=maxd=0;scaledv=scaledd=0
  for key,r in items.items():
   if key not in refq:continue
   ref=refq[key];v=complex_field(r,'value');d=complex_field(r,'derivative');rv=complex_field(ref,'value');rd=complex_field(ref,'derivative')
   if int(r['status'])==0 and int(ref['status'])==0:
    maxv=max(maxv,abs(v-rv)/max(abs(rv),qual['maximum']*1e-12));maxd=max(maxd,abs(d-rd)/max(abs(rd),qual['maximum']*1e-12))
    ri=int(r['row']);rowpeak=max(t[2] for t in F[ri]);slope=max(abs(F[ri][j+1][2]-F[ri][j][2])/(us[j+1]-us[j]) for j in range(m-1))
    scaledv=max(scaledv,abs(v-rv)/max(abs(rv),1e-3*rowpeak));scaledd=max(scaledd,abs(d-rd)/max(abs(rd),1e-3*slope))
  summary.append([b,n,len(items),outside,json.dumps(statuses),maxv,maxd,scaledv,scaledd])
 for row,(_,key) in selected.items():
  ref=refq[key];rv=complex_field(ref,'value');rd=complex_field(ref,'derivative')
  vals=[]
  for n in [10,20,30,40]:
   r=replays.get(('native',n),{}).get(key)
   vals.append(abs(complex_field(r,'value')-rv)/max(abs(rv),qual['maximum']*1e-12) if r else float('nan'))
  behavior.append([row,float(ref['u_re']),float(ref['u_im']),*vals])
  if row in rows:
   z=complex(float(ref['x_re']),float(ref['x_im']))
   for n in [10,20,30,40]:
    cc={int(r[3]):float(r[4]) for r in coeff if r[0]=='native' and r[1]==n and int(r[2])==row};t0=1+0j;t1=z
    for j in range(n+1):
     t=t0 if j==0 else t1 if j==1 else 2*z*t1-t0
     if j>=2:t0,t1=t1,t
     terms.append([row,n,j,z.real,z.imag,abs(cc[j]),abs(cc[j]*t)])
write('query_summary.csv',['backend','degree','queries','outside_real_interval','statuses','max_value_difference_from_native40','max_derivative_difference_from_native40','row_scaled_value_difference','row_scaled_derivative_difference'],summary)
write('complex_behavior.csv',['row','u_re','u_im','native10','native20','native30','native40'],behavior)
write('complex_terms.csv',['row','degree','order','x_re','x_im','abs_coefficient','abs_contribution'],terms)
if behavior:
 with (out/'complex_behavior.dat').open('w') as f:
  for r in behavior:f.write(' '.join(map(str,r))+'\n')
 script="set term svg size 1000,600 dynamic;set output 'complex_sensitivity.svg';set xlabel 'perpendicular row';set ylabel '|F_N-F_40| / max(|F_40|, peak*1e-12)';set logscale y;plot for [c=4:6] 'complex_behavior.dat' u 1:c w l t sprintf('native N%d',10*(c-3))\n"
 (out/'complex.gp').write_text(script);subprocess.run(['gnuplot','complex.gp'],cwd=out,check=True)
# Branch plot preserves raw results; qualification is reported separately.
branches=[r for r in results if r[4]=='scan' and r[0].startswith('scan_') and r[9] in ['root','scan']]
if branches:
 with (out/'branches.dat').open('w') as f:
  for b in ['legacy','shared','svd_index_control','native']:
   for n in [10,20,30,40]:
    for sel in [0,1]:
     rr=sorted([r for r in branches if r[1]==b and r[2]==n and r[3]==sel],key=lambda r:float(r[11]))
     if not rr:continue
     for r in rr:f.write(f'{r[11]} {r[12]} {r[13]}\n')
     f.write('\n\n')
 script="set term svg size 1100,550 dynamic;set output 'branches.svg';set multiplot layout 1,2;set xlabel 'k_parallel';set ylabel 'Re omega';plot 'branches.dat' u 1:2 w lp notitle;set ylabel 'Im omega';plot 'branches.dat' u 1:3 w lp notitle;unset multiplot\n"
 (out/'branches.gp').write_text(script);subprocess.run(['gnuplot','branches.gp'],cwd=out,check=True)
# Root/branch qualification and phase-independent polarization diagnostics.
root_records=[r for r in results if r[9] in ['root','scan']]
reference={}
for r in root_records:
 if r[1]=='legacy' and r[2]==30 and r[7]==4:
  reference[r[3],round(float(r[11]),12)]=r
quality=[]
for r in root_records:
 e=[complex(float(r[18+i]),float(r[21+i])) for i in range(3)] # prefix(9) + result fields E_re(9), E_im(12)
 norm=sum(abs(z)**2 for z in e);longitudinal=abs(e[2])**2/norm
 helicity=2*(e[0]*e[1].conjugate()).imag/max(abs(e[0])**2+abs(e[1])**2,1e-300)
 ref=reference.get((r[3],round(float(r[11]),12)))
 overlap='NA';dw=dg='NA'
 if ref:
  er=[complex(float(ref[18+i]),float(ref[21+i])) for i in range(3)]
  overlap=abs(sum(a.conjugate()*b for a,b in zip(e,er)))/math.sqrt(norm*sum(abs(z)**2 for z in er))
  dw=float(r[12])-float(ref[12]);dg=float(r[13])-float(ref[13])
 algebraic=int(r[14])<49 and float(r[17])<=1e-7
 state='matched numerical branch' if algebraic and isinstance(overlap,float) and overlap>=.95 else 'unqualified or no reference'
 quality.append([r[0],r[1],r[2],r[3],r[4],r[7],r[10],r[11],r[12],r[13],r[14],r[17],longitudinal,helicity,overlap,dw,dg,state])
write('root_qualification.csv',['case','backend','degree','selector','phase','ranks','point','kpar','omega_re','omega_im','iterations','scaled_residual','longitudinal_E_fraction','xy_helicity','overlap_legacy30','delta_omega_re','delta_gamma','qualification'],quality)
# Make the scientifically relevant diagonal imaginary response visible alongside global tensor errors.
rows33=[]
for r in tensors:
 if r[4]=='fixed' and r[7]==4 and r[5]==1 and r[6]==.075 and r[14]=='1' and r[15]=='3' and r[16]=='3':
  rows33.append([r[1],r[2],r[3],r[17],r[18]])
write('proton_chi33_seed.csv',['backend','degree','selector','real','imaginary'],rows33)
print('Root/branch records:',len(quality),flush=True)
# Explicit total susceptibilities, determinant deltas, and labelled dispersion panels.
rdata=list(csv.DictReader((out/'results.csv').open()))
tdata=list(csv.DictReader((out/'tensors.csv').open()))
refs={(r['degree'],r['selector'],r['probe'],r['initial_kpar'],r['ranks']):r for r in rdata if r['phase']=='fixed' and r['backend']=='legacy'}
deltas=[]
for r in rdata:
 k=tuple(r[v] for v in ['degree','selector','probe','initial_kpar','ranks'])
 if r['phase']=='fixed' and k in refs:
  ref=refs[k];d=complex(float(r['det_re']),float(r['det_im']))-complex(float(ref['det_re']),float(ref['det_im']))
  deltas.append([r['case'],r['backend'],d.real,d.imag,abs(d)])
write('determinant_differences.csv',['case','backend','delta_re','delta_im','absolute_complex_difference'],deltas)
total={}
for r in tdata:
 if r['species']=='0':continue
 k=(r['case'],r['label'],r['point'],r['i'],r['j']);total[k]=total.get(k,0j)+complex(float(r['value_re']),float(r['value_im']))
write('total_susceptibility.csv',['case','label','point','i','j','real','imaginary'],[list(k)+[v.real,v.imag] for k,v in total.items()])
groups=[('legacy',30),('shared',30),('svd_index_control',30),('native',10),('native',20),('native',30),('native',40)]
gp="set term svg size 1200,800 dynamic;set output 'branches.svg';set multiplot layout 2,2;set key font ',8';set xlabel 'k_parallel'\n"
for sel in ['0','1']:
 paths=[]
 for b,n in groups:
  rr=[r for r in rdata if r['phase']=='scan' and r['backend']==b and r['degree']==str(n) and r['selector']==sel and r['label'] in ['root','scan']]
  if not rr:continue
  fn=f'branch_{b}_{n}_{sel}.dat'
  with(out/fn).open('w') as f:
   for r in sorted(rr,key=lambda r:float(r['kpar'])):f.write(f"{r['kpar']} {r['omega_re']} {r['omega_im']}\n")
  paths.append((fn,f'{b} N{n}'))
 for col,title in [(2,'Re omega'),(3,'Im omega')]:
  gp+=f"set title 'selector {sel}';set ylabel '{title}';plot "+','.join(f"'{fn}' u 1:{col} w lp t '{name}'" for fn,name in paths)+'\n'
gp+='unset multiplot\n';(out/'branches.gp').write_text(gp)
if any(r['phase']=='scan' for r in rdata):subprocess.run(['gnuplot','branches.gp'],cwd=out,check=True)
# Singleton test assembly versus supported four-rank assembly, with raw deltas.
mpi_reference={}
for r in tdata:
 if r['ranks']=='4' and r['phase']=='fixed' and r['degree']=='30':
  mpi_reference[r['backend'],r['selector'],r['probe'],r['initial_kpar'],r['species'],r['i'],r['j']]=r
mpi_rows=[]
for r in tdata:
 key=tuple(r[k] for k in ['backend','selector','probe','initial_kpar','species','i','j'])
 if r['ranks']=='1' and r['phase']=='fixed' and key in mpi_reference:
  ref=mpi_reference[key];dr=float(r['value_re'])-float(ref['value_re']);di=float(r['value_im'])-float(ref['value_im'])
  mpi_rows.append([*key,r['value_re'],r['value_im'],ref['value_re'],ref['value_im'],dr,di])
write('mpi_tensor_differences.csv',['backend','selector','probe','kpar','species','i','j','singleton_re','singleton_im','four_rank_re','four_rank_im','delta_re','delta_im'],mpi_rows)

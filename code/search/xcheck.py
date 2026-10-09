# Independent-pipeline check: per-sector depth at fixed ephemeris in QLP (and TESS-SPOC) separately, plus per-sector blind TLS
import sys, os, subprocess, json, warnings, numpy as np; warnings.filterwarnings('ignore')
import lightkurve as lk
tic=int(sys.argv[1]); P=float(sys.argv[2]); T0=float(sys.argv[3]); D=float(sys.argv[4])/24
ts=f'{tic:016d}'; p=f'{ts[0:4]}/{ts[4:8]}/{ts[8:12]}/{ts[12:16]}'
sr=lk.search_lightcurve(f'TIC {tic}')
secs=sorted(set(int(m.split()[-1]) for m in sr.mission))
out={}
for prod in ('QLP','TESS-SPOC'):
  T=[];F=[];S=[]
  for s in secs:
    if prod=='QLP': url=f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{p}/hlsp_qlp_tess_ffi_s{s:04d}-{ts}_tess_v01_llc.fits'; fn=f'lc/{tic}/qlp_s{s}.fits'; cols=('kspsap_flux','det_flux')
    else: url=f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{p}/hlsp_tess-spoc_tess_phot_{ts}-s{s:04d}_tess_v1_lc.fits'; fn=f'lc/{tic}/spoc_s{s}.fits'; cols=('pdcsap_flux',)
    if not os.path.exists(fn): subprocess.run(['curl','-sf','-m','60','-o',fn,url])
    if not os.path.exists(fn) or os.path.getsize(fn)==0: continue
    lc=None
    for c in cols:
      try: lc=lk.read(fn,flux_column=c,quality_bitmask='hard'); break
      except Exception: pass
    if lc is None: continue
    lc=lc.remove_nans().normalize(); cad=np.nanmedian(np.diff(lc.time.value))
    t=lc.time.value; ph=((t-T0+0.5*P)%P)-0.5*P
    mask=np.abs(ph)<D
    lc=lc.flatten(window_length=int(0.75/cad)//2*2+1,mask=mask).remove_outliers(sigma_upper=3,sigma_lower=20)
    t=lc.time.value; f=np.asarray(lc.flux.value,float); ph=((t-T0+0.5*P)%P)-0.5*P
    it=np.abs(ph)<D/2*0.8; oo=(np.abs(ph)>D)&(np.abs(ph)<3*D)
    # empirical error: same statistic at 300 random phases
    rs=[]
    rng=np.random.default_rng(1)
    for x in rng.uniform(0.15,0.85,300):
      p2=((t-T0-x*P+0.5*P)%P)-0.5*P; a=np.abs(p2)<D/2*0.8; b=(np.abs(p2)>D)&(np.abs(p2)<3*D)
      if a.sum()>2: rs.append(np.median(f[b])-np.mean(f[a]))
    d=np.median(f[oo])-np.mean(f[it])
    out[f'{prod}_S{s}']=dict(cad_min=round(cad*1440,1),depth_ppm=round(d*1e6),err_emp_ppm=round(np.std(rs)*1e6),n_in=int(it.sum()))
    T.append(t);F.append(f);S.append(np.full(len(t),s))
  if T:
    t=np.concatenate(T);f=np.concatenate(F)
    np.savez(f'lc/{tic}_{prod}_masked.npz',t=t,f=f,s=np.concatenate(S))
    ph=((t-T0+0.5*P)%P)-0.5*P; it=np.abs(ph)<D/2*0.8; oo=(np.abs(ph)>D)&(np.abs(ph)<3*D)
    ds=[v['depth_ppm'] for k,v in out.items() if k.startswith(prod)]; es=[v['err_emp_ppm'] for k,v in out.items() if k.startswith(prod)]
    w=1/np.array(es,float)**2; m=np.sum(np.array(ds)*w)/w.sum(); e=1/np.sqrt(w.sum())
    out[f'{prod}_combined']=dict(depth_ppm=round(m),err_ppm=round(e),snr=round(m/e,2),chi2=round(float(np.sum((np.array(ds)-m)**2*w)),1),dof=len(ds)-1)
json.dump(out,open(f'xcheck_{tic}.json','w'),indent=1)
for k,v in out.items(): print(k,v)

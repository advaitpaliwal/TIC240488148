# Re-run TRICERATOPS NFPP/FPP with ruled-out neighbours removed (target.remove_star). Same data and settings as hunt-ffi/trc/run_fpp.py.
# usage: trc_rerun.py TAG SEED ID1,ID2,...   (empty list = baseline)
import warnings; warnings.filterwarnings('ignore')
import sys, os, json, pickle, time, numpy as np
H='work/hunt-ffi/'
sys.path.insert(0,H+'trc')
def main():
    import gaia_patch
    import triceratops.triceratops as tr
    tag,k=sys.argv[1],int(sys.argv[2]); drop=[x for x in sys.argv[3].split(',') if x] if len(sys.argv)>3 else []
    tic,P,T0=240488148,1.512529,1712.925
    _p=H+f'lc/{tic}_TESS-SPOC_masked.npz'; z=np.load(_p if os.path.exists(_p) else H+f'lc/{tic}_QLP_masked.npz'); t,f=z['t'],z['f']
    tg=pickle.load(open(H+f'trc/out/target_{tic}.pkl','rb'))
    keep=~tg.stars['ID'].isin([str(d) for d in drop]).values
    tg.pix_coords=[pc[keep] for pc in tg.pix_coords]   # keep per-sector pixel coords aligned with .stars
    tg.remove_star(drop); tg.stars=tg.stars.reset_index(drop=True)
    ph=t-T0-np.round((t-T0)/P)*P; w=np.abs(ph)<0.25
    edges=np.arange(-0.25,0.25+1e-9,10/1440); idx=np.digitize(ph[w],edges)
    tb=np.array([ph[w][idx==i].mean() for i in np.unique(idx)]); fb=np.array([f[w][idx==i].mean() for i in np.unique(idx)])
    dur=0.0625; oot=np.abs(tb)>dur; err=float(np.std(fb[oot])); depth=float(np.median(fb[oot])-fb[np.abs(tb)<dur/2*0.8].mean())
    tg.calc_depths(tdepth=depth,all_ap_pixels=tg._aps)
    np.random.seed(1000+k); t0=time.time()
    tg.calc_probs(time=tb,flux_0=fb,flux_err_0=err,P_orb=P,contrast_curve_file=None,filt='TESS',N=1000000,parallel=True,exptime=30/1440,nsamples=20,verbose=0)
    tg.probs.to_csv(f'out/trc_probs_{tag}_{k}.csv',index=False)
    json.dump(dict(tag=tag,run=k,dropped=drop,n_stars=len(tg.stars),depth_ppm=depth*1e6,FPP=float(tg.FPP),NFPP=float(tg.NFPP),sec=time.time()-t0),open(f'out/trc_{tag}_{k}.json','w'))
    print(tag,k,len(drop),tg.FPP,tg.NFPP,round(time.time()-t0),flush=True)
if __name__=='__main__': main()

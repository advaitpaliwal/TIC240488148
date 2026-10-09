# TRICERATOPS FPP for TIC 240488148, same LC/binning/depth procedure as ../../hunt-ffi/trc/run_fpp.py and the same target pickle.
# usage: run_fpp.py MODE RUN_INDEX ; MODE in {none, cc50, cc90, molusc}
import sys, json, pickle, os, time, numpy as np, warnings; warnings.filterwarnings('ignore')
import gaia_patch
import triceratops.triceratops as tr
import types, pandas as pd, triceratops.marginal_likelihoods as ml
class _F:  # pandas>=3 copy-on-write makes Series.values read-only; TRICERATOPS 1.1.0 writes into it (molusc branch only)
    def __init__(s, d): s.d = d
    def __getitem__(s, k):
        if isinstance(k, str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
        return _F(s.d[np.asarray(k)])
_rc = ml.read_csv
ml.read_csv = lambda f, *a, **kw: _F(_rc(f, *a, **kw)) if 'molusc' in str(f) else _rc(f, *a, **kw)
mode, k = sys.argv[1], int(sys.argv[2])
tic, P, T0, dur = 240488148, 1.512529, 1712.925, 0.06
z = np.load(os.path.expanduser(f'work/hunt-ffi/lc/{tic}_TESS-SPOC_masked.npz')); t, f = z['t'], z['f']
tg = pickle.load(open('out/target_240488148_from_hunt-ffi.pkl', 'rb'))
ph = t - T0 - np.round((t - T0) / P) * P
w = np.abs(ph) < 0.25
edges = np.arange(-0.25, 0.25 + 1e-9, 10 / 1440)
idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)])
fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
oot = np.abs(tb) > dur
err = float(np.std(fb[oot]))
depth = float(np.median(fb[oot]) - fb[np.abs(tb) < dur / 2 * 0.8].mean())
tg.calc_depths(tdepth=depth, all_ap_pixels=tg._aps)
cc = {'none': None, 'cc50': '../cc_gaia_dr2completeness_50pct_G.csv', 'cc90': '../cc_gaia_dr2completeness_90pct_G.csv', 'molusc': None}[mode]
mf = '../molusc_out/molusc_for_trc_N1e6.csv' if mode == 'molusc' else None
np.random.seed(1000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=cc, filt='Vis',
              N=1000000, parallel=True, exptime=30 / 1440, nsamples=20, verbose=0, molusc_file=mf)
tag = f'{mode}_{k}'
tg.probs.to_csv(f'out/probs_{tag}.csv', index=False)
json.dump(dict(mode=mode, run=k, depth_ppm=depth * 1e6, err_ppm=err * 1e6, nbins=len(tb), cc=cc, molusc=mf,
               FPP=float(tg.FPP), NFPP=float(tg.NFPP), sec=time.time() - t0), open(f'out/res_{tag}.json', 'w'))
print(tag, tg.FPP, tg.NFPP, round(time.time() - t0), flush=True)

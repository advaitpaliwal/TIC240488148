# Adapted from ../../triceratops/{setup_target,run_fpp}.py for FFI-only stars.
# usage: run_fpp.py TIC P T0 RUN_INDEX
import sys, json, pickle, os, time, numpy as np, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
import gaia_patch
import triceratops.triceratops as tr
tic, P, T0, k = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
_p = f'../lc/{tic}_TESS-SPOC_masked.npz'; z = np.load(_p if os.path.exists(_p) else f'../lc/{tic}_QLP_masked.npz'); t, f, sid = z['t'], z['f'], z['s']
secs = np.unique(sid).astype(int)
pk = f'out/target_{tic}.pkl'
if not os.path.exists(pk):
    print('building target',flush=True); tg = tr.target(ID=tic, sectors=secs); print('target ok',flush=True)
    aps = []
    for i, s in enumerate(secs):
        fn = f'../lc/{tic}/spoc_s{s}.fits'
        ap = None
        if os.path.exists(fn):
            with fits.open(fn) as h:
                a = h['APERTURE'].data.astype(int); c0 = h['APERTURE'].header['CRVAL1P']; r0 = h['APERTURE'].header['CRVAL2P']
            rows, cols = np.nonzero(np.bitwise_and(a, 2))
            if len(rows): ap = np.column_stack([cols + c0, rows + r0])
        if ap is None:  # QLP: no aperture in file; 5x5 centred on target (TRICERATOPS default)
            tp = np.round(tg.pix_coords[i][0])
            ap = np.array([np.repeat(np.arange(tp[0]-2, tp[0]+3), 5), np.tile(np.arange(tp[1]-2, tp[1]+3), 5)]).T
        aps.append(ap)
    tg._aps = aps
    pickle.dump(tg, open(pk, 'wb'))
tg = pickle.load(open(pk, 'rb'))
ph = t - T0 - np.round((t - T0) / P) * P
w = np.abs(ph) < 0.25
edges = np.arange(-0.25, 0.25 + 1e-9, 10 / 1440)
idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)])
fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
dur = float(sys.argv[5]) if len(sys.argv) > 5 else 0.06
oot = np.abs(tb) > dur
err = float(np.std(fb[oot]))
depth = float(np.median(fb[oot]) - fb[np.abs(tb) < dur / 2 * 0.8].mean())
tg.calc_depths(tdepth=depth, all_ap_pixels=tg._aps)
np.random.seed(1000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=None, filt='TESS',
              N=1000000, parallel=True, exptime=30 / 1440, nsamples=20, verbose=0)
tag = f'{tic}_{k}'
tg.probs.to_csv(f'out/probs_{tag}.csv', index=False)
json.dump(dict(tic=tic, run=k, P=P, T0=T0, depth_ppm=depth * 1e6, err_ppm=err * 1e6, nbins=len(tb),
               FPP=float(tg.FPP), NFPP=float(tg.NFPP), sec=time.time() - t0), open(f'out/res_{tag}.json', 'w'))
print(tag, tg.FPP, tg.NFPP, round(time.time() - t0))

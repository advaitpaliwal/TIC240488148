# Glue: TESS-SPOC PDCSAP -> wotan detrend (transit masked) -> batman + emcee fit for TIC 240488148.01
import numpy as np, emcee, batman, json, sys
from astropy.io import fits
from wotan import flatten, transit_mask
from scipy.optimize import minimize
np.random.seed(42)
import os
LC = os.path.join(os.environ.get('TIC_LC_DIR', 'lc'), 'spoc_s%d.fits')  # TESS-SPOC FFI light curves from MAST
SECS = [15, 16, 56, 76, 82, 83]
P0, T00, DUR0 = 1.512529, 1712.9250, 1.5/24
W = 3*DUR0  # half-window kept around each transit for the fit

# Stellar priors, TIC v8.2 (MAST, _sources/tic240488148_mast.csv)
RHO_SUN = 1.40894  # g/cm3
RHO, E_RHO = 2.27506*RHO_SUN, 0.705757*RHO_SUN
# Claret 2017 TESS quadratic LD, ATLAS LSM, Z=0, xi=2; bilinear in Teff (4250,4500) and logg (4.5,5.0) at 4403 K, 4.622
def interp(v4250_45, v4500_45, v4250_50, v4500_50, T=4403., g=4.62235):
    x=(T-4250)/250; y=(g-4.5)/0.5
    return (1-x)*(1-y)*v4250_45 + x*(1-y)*v4500_45 + (1-x)*y*v4250_50 + x*y*v4500_50
import csv
rows=[r for r in csv.DictReader(open('_sources/claret2017_table25_sel.csv')) if float(r['xi'])==2.0]
g={(float(r['Teff']),float(r['logg'])):(float(r['aLSM']),float(r['bLSM'])) for r in rows}
U1=interp(g[(4250,4.5)][0],g[(4500,4.5)][0],g[(4250,5.0)][0],g[(4500,5.0)][0])
U2=interp(g[(4250,4.5)][1],g[(4500,4.5)][1],g[(4250,5.0)][1],g[(4500,5.0)][1])
E_U=0.1

T, F, E, S = [], [], [], []
for s in SECS:
    d = fits.open(LC % s)[1].data
    m = (d['QUALITY'] == 0) & np.isfinite(d['PDCSAP_FLUX']) & np.isfinite(d['TIME'])
    t, f, e = d['TIME'][m], d['PDCSAP_FLUX'][m], d['PDCSAP_FLUX_ERR'][m]
    e = e/np.nanmedian(f); f = f/np.nanmedian(f)
    msk = transit_mask(time=t, period=P0, duration=2.5*DUR0, T0=T00)
    trd = flatten(t[~msk], f[~msk], method='biweight', window_length=0.75, return_trend=True)[1]
    good = np.isfinite(trd); fl = np.interp(t, t[~msk][good], trd[good])
    # drop points farther than 0.25 d from any out-of-transit point used by the trend (gap edges)
    near = np.abs(t[:, None] - t[~msk][good][None, ::max(1, good.sum()//4000)]).min(1) < 0.25
    t, f, e = t[near], f[near]/fl[near], e[near]/fl[near]
    oot = ~transit_mask(time=t, period=P0, duration=2.5*DUR0, T0=T00)
    r = f - 1; sig = 1.4826*np.median(np.abs(r[oot]-np.median(r[oot])))
    keep = (r < 5*sig) & (r > -10*sig)  # clip upward outliers / extreme dips only
    ph = ((t - T00 + 0.5*P0) % P0) - 0.5*P0
    w = keep & (np.abs(ph) < W)
    T.append(t[w]); F.append(f[w]); E.append(e[w]); S.append(np.full(w.sum(), s))
    print(s, w.sum(), 'oot rms ppm', round(sig*1e6))
T, F, E, S = map(np.concatenate, (T, F, E, S))
long = np.isin(S, [15, 16])
np.savez('fitdata.npz', t=T, f=F, e=E, s=S)

TREF_N = int(np.round((np.median(T) - T00)/P0))  # reference epoch near middle of baseline
TREF0 = T00 + TREF_N*P0
p = batman.TransitParams(); p.t0=0.; p.per=P0; p.rp=0.02; p.a=8.; p.inc=88.; p.ecc=0.; p.w=90.; p.limb_dark='quadratic'; p.u=[U1,U2]
mL = batman.TransitModel(p, T[long], supersample_factor=15, exp_time=1800/86400.)
mS = batman.TransitModel(p, T[~long], supersample_factor=3, exp_time=200/86400.)

def a_rs(rho, P):  # rho g/cm3, P days -> a/R*
    G=6.674e-8; return ((G*rho*1e0*(P*86400)**2)/(3*np.pi))**(1/3)

names = ['T0', 'P', 'k', 'b', 'rho', 'q1', 'q2', 'f0', 'lnsL', 'lnsS']
def unpack(th):
    T0, P, k, b, rho, q1, q2, f0, lsL, lsS = th
    a = a_rs(rho, P); inc = np.degrees(np.arccos(b/a))
    u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1-2*q2)
    return T0, P, k, b, rho, a, inc, u1, u2, f0, lsL, lsS
def model(th, which=None):
    T0, P, k, b, rho, a, inc, u1, u2, f0, _, _ = unpack(th)
    p.t0=T0; p.per=P; p.rp=k; p.a=a; p.inc=inc; p.u=[u1,u2]
    out = np.empty_like(F); out[long]=mL.light_curve(p); out[~long]=mS.light_curve(p)
    return out*f0
def lnprob(th):
    T0, P, k, b, rho, q1, q2, f0, lsL, lsS = th
    if not (0 < k < 0.1 and 0 <= b < 1+k and 0.1 < rho < 20 and 0 < q1 < 1 and 0 < q2 < 1 and -12 < lsL < -4 and -12 < lsS < -4): return -np.inf
    if abs(T0 - TREF0) > 0.05 or abs(P - P0) > 0.001: return -np.inf
    _, _, _, _, _, a, inc, u1, u2, _, _, _ = unpack(th)
    if b >= a: return -np.inf
    lp = -0.5*((rho-RHO)/E_RHO)**2 - 0.5*((u1-U1)/E_U)**2 - 0.5*((u2-U2)/E_U)**2
    r = F - model(th)
    var = E**2 + np.where(long, np.exp(2*lsL), np.exp(2*lsS))
    return lp - 0.5*np.sum(r**2/var + np.log(var))

if __name__ == '__main__':
    q1_0=(U1+U2)**2; q2_0=U1/(2*(U1+U2))
    x0 = np.array([TREF0, P0, 0.0193, 0.3, RHO, q1_0, q2_0, 1.0, -9, -9])
    nll = lambda x: -lnprob(x)
    best = x0
    for _ in range(3):
        r = minimize(nll, best, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-3)); best = r.x
    print('MAP', dict(zip(names, best)), -r.fun)
    nw, nd = 40, len(best)
    sc = np.array([2e-4, 2e-6, 1e-3, 0.05, 0.2, 0.02, 0.02, 1e-5, 0.05, 0.05])
    pos = []
    while len(pos) < nw:
        c = best + sc*np.random.randn(nd)
        if np.isfinite(lnprob(c)): pos.append(c)
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 15000
    be = emcee.backends.HDFBackend('chain.h5') if False else None
    sam = emcee.EnsembleSampler(nw, nd, lnprob)
    sam.run_mcmc(np.array(pos), nsteps, progress=False)
    try: tau = sam.get_autocorr_time(tol=0)
    except Exception as ex: tau = np.full(nd, np.nan)
    print('tau', dict(zip(names, np.round(tau, 1))))
    burn = nsteps//3
    ch = sam.get_chain(discard=burn, thin=10, flat=True)
    np.savez('chain.npz', chain=ch, tau=tau, names=names, nsteps=nsteps, burn=burn, nw=nw,
             acc=sam.acceptance_fraction, TREF_N=TREF_N, U1=U1, U2=U2, best=best)
    print('acc', np.mean(sam.acceptance_fraction), 'n', len(ch))

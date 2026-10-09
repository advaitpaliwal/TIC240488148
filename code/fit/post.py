# Post-process chain.npz -> derived parameters, corner plot, folded model plot, RESULT numbers
import numpy as np, sys, json, batman, corner
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os; sys.argv = ['x']; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import fit as Fm
np.random.seed(1)
c = np.load('chain.npz', allow_pickle=True); ch = c['chain']; N = len(ch)
T0, P, k, b, rho, q1, q2, f0, lsL, lsS = ch.T
a = Fm.a_rs(rho, P); inc = np.degrees(np.arccos(b/a))
u1 = 2*np.sqrt(q1)*q2; u2 = np.sqrt(q1)*(1-2*q2)
# total duration T14 (circular)
T14 = P/np.pi*np.arcsin(np.sqrt((1+k)**2-b**2)/(a*np.sin(np.radians(inc))))*24
# depth: max model depth of the unbinned model
def depth(i):
    p = batman.TransitParams(); p.t0=0; p.per=P[i]; p.rp=k[i]; p.a=a[i]; p.inc=inc[i]; p.ecc=0; p.w=90; p.limb_dark='quadratic'; p.u=[u1[i],u2[i]]
    return 1-batman.TransitModel(p, np.array([0.])).light_curve(p)[0]
dep = np.array([depth(i) for i in range(N)])*1e6
# stellar: independent Gaussian draws of TIC R*, M*, Teff, L (TIC v8.2)
Rs = np.random.normal(0.67187, 0.0571432, N); Ms = np.random.normal(0.69, 0.0827815, N)
Teff = np.random.normal(4403., 127.12, N); L = np.random.normal(0.152848229, 0.008660045, N)
Rp = k*Rs*109.076  # R_sun / R_earth (IAU nominal)
sma_au = (Ms*(P/365.25)**2)**(1/3)                # Kepler III with TIC mass
sma_au_rho = a*Rs*0.00465047                      # from fitted a/R* x TIC R*
Teq = Teff*np.sqrt(1/(2*a))                       # zero albedo, full redistribution
S = L/sma_au**2                                   # insolation, S_earth, from TIC luminosity
BJD = 2457000.
NREF = int(c['TREF_N'])
# ephemeris prediction for Nov 2026 (BJD 2461345.5 = 2026-11-01 0h UTC approx; 2461375.5 = 2026-12-01)
def nidx(bjd): return np.round((bjd-BJD-np.median(T0))/np.median(P))
nn = np.arange(nidx(2461345.5), nidx(2461375.5)+1)
pred = np.array([T0 + n*P for n in nn])
pred_sig_min = pred.std(1)*1440
def q(x): m, lo, hi = np.percentile(x, [50, 15.865, 84.135]); return m, hi-m, m-lo
out = {}
for nm, x in [('P_d', P), ('T0_BJD', T0+BJD), ('depth_ppm', dep), ('T14_h', T14), ('k', k), ('a_Rs', a), ('b', b), ('inc_deg', inc),
              ('rho_star_gcc', rho), ('u1', u1), ('u2', u2), ('Rp_Re', Rp), ('Teq_K', Teq), ('S_Se', S), ('a_au', sma_au), ('a_au_from_aRs', sma_au_rho),
              ('jitter30_ppm', np.exp(lsL)*1e6), ('jitter200_ppm', np.exp(lsS)*1e6)]:
    out[nm] = q(x)
out['b_95pct_upper'] = float(np.percentile(b, 95))
out['N_samples'] = N; out['tau'] = dict(zip(c['names'].tolist(), np.round(c['tau'], 1).tolist())); out['nsteps'] = int(c['nsteps']); out['burn'] = int(c['burn'])
out['acc_mean'] = float(np.mean(c['acc'])); out['epoch_ref_n_from_search_T0'] = NREF
out['pred_nov2026'] = dict(first_epoch=int(nn[0]), first_BJD=float(np.median(pred[0])+BJD), sig_min_range=[float(pred_sig_min.min()), float(pred_sig_min.max())])
# P-T0 correlation
out['corr_T0_P'] = float(np.corrcoef(T0, P)[0, 1])
json.dump(out, open('result.json', 'w'), indent=1, default=float)
for kx, v in out.items(): print(kx, v)

# corner plot
lab = [r'$T_0-%.4f$ [min]' % (np.median(T0)), r'$P-%.6f$ [s]' % np.median(P), r'$R_p/R_\star$', r'$b$', r'$\rho_\star$ [g cm$^{-3}$]', r'$u_1$', r'$u_2$']
X = np.vstack([(T0-np.median(T0))*1440, (P-np.median(P))*86400, k, b, rho, u1, u2]).T
fig = corner.corner(X, labels=lab, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt='.3f', title_kwargs=dict(fontsize=9))
fig.suptitle('TIC 240488148.01  batman + emcee  (TESS-SPOC S15,16,56,76,82,83)', y=1.01)
fig.savefig('corner_TIC240488148.png', dpi=110, bbox_inches='tight')

# folded model-over-data plot
med = np.median(ch, 0)
ph = ((Fm.T - np.median(T0) + 0.5*np.median(P)) % np.median(P) - 0.5*np.median(P))*24
fn = Fm.F/np.median(f0)
fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
tt = np.linspace(-0.2, 0.2, 2000)/1  # days
p = batman.TransitParams(); p.t0=0; p.per=np.median(P); p.rp=np.median(k); p.a=np.median(a); p.inc=np.median(inc); p.ecc=0; p.w=90; p.limb_dark='quadratic'; p.u=[np.median(u1), np.median(u2)]
mod_s = batman.TransitModel(p, tt, supersample_factor=3, exp_time=200/86400).light_curve(p)
mod_l = batman.TransitModel(p, tt, supersample_factor=15, exp_time=1800/86400).light_curve(p)
for i in np.random.choice(N, 50):
    pp = batman.TransitParams(); pp.t0=0; pp.per=P[i]; pp.rp=k[i]; pp.a=a[i]; pp.inc=inc[i]; pp.ecc=0; pp.w=90; pp.limb_dark='quadratic'; pp.u=[u1[i], u2[i]]
    ax[0].plot(tt*24, (batman.TransitModel(pp, tt, supersample_factor=3, exp_time=200/86400).light_curve(pp)-1)*1e6, c='C1', alpha=0.06, lw=1)
ax[0].plot(Fm.T[:0], Fm.F[:0])
ax[0].plot(ph, (fn-1)*1e6, '.', c='0.82', ms=1.5, zorder=0, label='TESS-SPOC (30 min and 200 s)')
bins = np.arange(-4.5, 4.51, 0.25); bc = (bins[1:]+bins[:-1])/2; idx = np.digitize(ph, bins)
mb = np.array([np.mean(fn[idx == j]) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))])
eb = np.array([np.std(fn[idx == j])/np.sqrt((idx == j).sum()) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))])
ax[0].errorbar(bc, (mb-1)*1e6, eb*1e6, fmt='o', c='k', ms=4, label='15-min bins')
ax[0].plot(tt*24, (mod_s-1)*1e6, c='C1', lw=2, label='median model (200 s)')
ax[0].plot(tt*24, (mod_l-1)*1e6, c='C0', lw=1.2, ls='--', label='median model (30 min smeared)')
ax[0].set_ylim(-900, 500); ax[0].set_ylabel('relative flux [ppm]'); ax[0].legend(fontsize=7, loc='lower left')
ax[0].set_title('TIC 240488148.01  P = %.6f d   Rp = %.2f R$_\\oplus$' % (np.median(P), np.median(Rp)))
# residual: data minus per-point model
res = Fm.F - Fm.model(med)
rb = np.array([np.mean(res[idx == j]) if (idx == j).sum() > 3 else np.nan for j in range(1, len(bins))])
ax[1].errorbar(bc, rb*1e6, eb*1e6, fmt='o', c='k', ms=3); ax[1].axhline(0, c='C1')
ax[1].set_ylim(-300, 300); ax[1].set_ylabel('resid. [ppm]'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_xlim(-4.5, 4.5)
plt.tight_layout(); plt.savefig('folded_model_TIC240488148.png', dpi=130)

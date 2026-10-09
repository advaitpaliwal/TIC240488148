# TIC 240488148.01: final transit fit (2026-10-08)

Method (verified, code in this folder): TESS-SPOC PDCSAP FFI light curves for S15, S16, S56, S76, S82, S83 (QUALITY==0; PDCSAP already corrects for crowding, CROWDSAP 0.88-0.94 from headers). Detrended with wotan biweight (0.75 d window) fitted to out-of-transit points only, with transits masked at 2.5x duration. Only ±4.5 h windows around each transit are kept (10,263 points). The model is batman, quadratic limb darkening, circular orbit. It is supersampled 15x over 1800 s for S15-16 and 3x over 200 s for the rest. Sampling is emcee, 40 walkers x 8000 steps with 1/3 burn-in. The maximum autocorrelation time is 237 steps, so the run is about 34τ; mean acceptance 0.33; 21,320 thinned samples. Fitted parameters are T0, P, Rp/R*, b, ρ*, Kipping q1 and q2, flux offset, and one jitter term for each cadence.
Priors: ρ* ~ N(3.21, 0.99) g/cm³ (TIC v8.2 ρ = 2.275 ± 0.706 solar). u1, u2 ~ N(0.494, 0.1) and N(0.171, 0.1), interpolated from the Claret (2017) TESS ATLAS tables at Teff 4403 K, logg 4.62, [M/H]=0. Inputs are saved in `_sources/`.
Stellar inputs (verified, MAST TIC v8.2, `_sources/tic240488148_mast.csv`): R* = 0.672 ± 0.057 R☉, M* = 0.690 ± 0.083 M☉, Teff = 4403 ± 127 K, L = 0.1528 ± 0.0087 L☉, logg = 4.622 ± 0.103.

## Posterior (median, +/- 68% interval)

| Parameter | Value |
|---|---|
| Period P | 1.5125149 +0.0000030 / -0.0000027 d (±0.25 s) |
| Mid-transit T0 | 2460393.34348 +0.00153 / -0.00115 BJD_TDB (±1.9 min; epoch chosen near the middle of the data, so the P-T0 correlation is only 0.54) |
| Depth (maximum model depth) | 411 +28 / -28 ppm |
| Duration T14 | 1.41 +0.07 / -0.06 h |
| Rp/R* | 0.01873 +0.00089 / -0.00085 |
| a/R* | 7.59 +0.54 / -0.76 |
| Impact parameter b | 0.42 +0.18 / -0.24 (95% upper limit 0.69) |
| Inclination | 86.8 +1.8 / -1.8 deg |
| ρ* (posterior) | 3.62 +0.83 / -0.98 g/cm³ |
| u1, u2 (posterior) | 0.50 ± 0.10, 0.17 ± 0.10 (prior-dominated) |
| Planet radius | 1.37 +0.14 / -0.13 R⊕ (includes the TIC R* error) |
| Equilibrium temperature | 1133 +68 / -51 K (A=0, full redistribution, Teff·sqrt(R*/2a)) |
| Insolation | 294 +32 / -27 S⊕ (TIC L / a²) |
| Semi-major axis | 0.02279 +0.00088 / -0.00095 AU (Kepler III with TIC M*); from a/R* × R*: 0.0235 ± 0.0029 AU |
| Jitter | 46 (+58/-34) ppm at 30 min, 22 (+31/-12) ppm at 200 s |

ExoFOP-style 1σ errors (symmetrised): P 1.5125149 ± 0.0000029 d; epoch 2460393.3435 ± 0.0013 (BJD_TDB); depth 411 ± 28 ppm; duration 1.41 ± 0.06 h; Rp/R* 0.0187 ± 0.0009; a/R* 7.6 ± 0.6; b 0.42 ± 0.21; inc 86.8 ± 1.8; Rp 1.37 ± 0.13 R⊕; Teq 1133 ± 60 K; S 294 ± 29 S⊕; a 0.0228 ± 0.0009 AU.

## Predicted transits, Nov 2026 (verified from the posterior, linear ephemeris)
1σ timing uncertainty is 4.0-4.1 min for every transit from 1 to 30 Nov 2026 (epochs 630-649 counted from T0 above). Examples: 2461346.2279 (2026-11-01 17:26 UTC), 2461347.7404 (11-03 05:44), 2461349.2529 (11-04 18:03), ..., 2461374.9656 (11-30 11:09), all ±4 min. The 3σ window is ±12 min.

## Notes
- Verified: the period is 1.4e-5 d shorter than the search value of 1.512529 d (about 5σ). With the search period, a per-sector timing scan (`ttv_check.py`) put S15/S16 transits 14 and 30 min late and S82/S83 on time. The fitted period removes that drift. Use the fitted ephemeris, not the search one.
- Verified: the 411 ppm fitted depth is deeper than the 374 ppm box-search depth because it is the limb-darkened model's central depth, not a mean in-transit level. (Rp/R*)² = 351 ppm.
- Inferred: the transit-derived ρ* (3.6 ± 0.9) agrees with TIC (3.2 ± 1.0), which is consistent with a planet on a near-circular orbit around the target.
- Inferred: b is mostly constrained by the ρ* prior, because b and ρ* are strongly degenerate (see corner plot).
- An earlier test run (superseded) had a detrending bug that removed in-transit points. The numbers above come from the fixed `fit.py`.

Files: `fit.py`, `post.py`, `ttv_check.py`, `chain.npz`, `result.json`, `run.log`, `corner_TIC240488148.png`, `folded_model_TIC240488148.png`, `_sources/`.

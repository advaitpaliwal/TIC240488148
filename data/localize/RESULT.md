# TIC 240488148 (P 1.512529 d, ~366 ppm): which star does the dip come from?

Run on 2026-10-08 in `~/Scratch/planet-hunt/tic240488148-localize/` with its own venv (`.venv`: tess-localize 0.7.6, triceratops 1.1.0, lightkurve 2.6.0, pandas 3.0.6, setuptools<81). Pixels are TESScut 15x15 FFI stamps for all six sectors (`cut/tesscut/`). Gaia DR3 comes from VizieR I/355/gaiadr3, with proper motion carried forward to each sector (`gaia_field.csv`, `common.py`).

## Answer

- **The signal does not come from TIC 240488141 or from any catalogued star 20" or more from the target.** In Sectors 15 and 16, each of those stars would need a flux drop of 6 to 236 e/s at its own position to produce the dip. The difference images show about 0 ± 0.2–0.7 e/s at every one of them. A forced fit at each star's required depth is worse than the target-only fit by Δχ² = 80 to 500,000 (two sectors summed). For TIC 240488141 specifically: Δχ² = 5,346; it would need a 20.9 / 10.2 e/s drop and shows 1.02 ± 0.34 / −0.05 ± 0.24 e/s. **[verified: `out/neighbour_tests.csv`, `post.py`]**
- **With those 26 stars removed, TRICERATOPS NFPP is 7.0e-9, 7.5e-9 and 6.9e-9 over three seeds (bar: NFPP < 0.001). That passes. FPP stays at 0.076–0.083, which fails the usual FPP < 0.015 bar.** The no-drop baseline rerun with the same script reproduces the original exactly: FPP 0.0868, NFPP 0.00376. **[verified: `out/trc_*.json`]**
- **Three faint stars within 12" cannot be excluded from pixels:** Gaia G 18.8–19.1 (TIC 2007804809 at 3.7", TIC 2007804808 at 7.3", TIC 240488136 at 11.2"). Each would need a 65–74% eclipse. They stayed in the TRICERATOPS rerun, which already gives them about 0 probability. **[verified]**
- **TESS-Localize does not cleanly place the source on the target.** In the 2-minute sectors (56, 76, 82, 83) its best-fit position lands 7–10" from the target, in among the target and those same three faint stars. In the two 30-minute sectors it lands 22–25" away at low or negative SNR, so those two fits are not usable. It never favours TIC 240488141 (relative likelihood 0 in all six sectors). **[verified numbers; my reading of why is inferred, see the caveats]**

## 1. TESS-Localize, per sector

Method: harmonics n/P for n = 1…25, 3 PCA components, 3x3 aperture on the target, Gaia stars down to G < 21 (`loc.py`, `out/loc_s*.json`, `out/loc_s*_starfit.csv`, heatmaps in `out/loc_s*.png`). RL means relative likelihood (normalised over 1,890 Gaia sources). Distances are in pixels (21"/px).

| Sector (cadence) | Fit offset from target | Fit error x/y (px) | Summed-height SNR | Top Gaia star (sep, G, RL) | Target RL (dist px) | TIC 240488141 RL (dist px) |
|---|---|---|---|---|---|---|
| 15 (30 min) | 22.3" | 0.11 / 0.16 | 3.8 | …1880448 (29.5", G19.5, 0.51) | 0.00017 (1.19) | 0 (2.86) |
| 16 (30 min) | 24.6" | 0.10 / 0.09 | **−4.7** | …9510272 (22.5", G18.0, 0.41) | 0 (1.26) | 0 (1.73) |
| 56 (2 min) | 9.2" | 0.10 / 0.09 | 8.8 | …6697216 (7.2", G18.8 = TIC 2007804808, 0.59) | 0.0029 (0.55) | 0 (2.21) |
| 76 (2 min) | 7.4" | 0.08 / 0.07 | 8.5 | …7022336 (11.3", G19.1 = TIC 240488136, 0.78) | **0.21** (0.26) | 0 (2.34) |
| 82 (2 min) | 8.4" | 0.08 / 0.10 | 7.9 | …6700288 (3.9", G18.9 = TIC 2007804809, 0.54) | 0.070 (0.42) | 0 (1.93) |
| 83 (2 min) | 9.7" | 0.10 / 0.12 | 7.9 | …6697216 (7.2", G18.8, 0.51) | 0.0018 (0.58) | 0 (2.18) |

- [verified] TIC 240488141 and every star of G < 16 other than the target get relative likelihood 0 in all six sectors.
- [inferred] TESS-Localize models a periodic signal as sinusoids. A 1.5-hour dip that fills about 4% of the orbit spreads its power across many low-amplitude harmonics, and the package's extrinsic error model was calibrated on binaries and pulsators. Its quoted errors (about 0.1 px, roughly 2") are therefore probably too small for a 366 ppm transit. The fact that all four 2-minute fits sit 0.26–0.58 px from the target but never on it looks like systematic bias, not an off-target source. The two 30-minute fits (SNR 3.8 and −4.7, the latter a negative summed amplitude) are noise.
- [verified] Because those errors are so tight, the "best" star in each sector is whichever faint G~19 star sits nearest the fit. Which one that is changes from sector to sector (TIC 2007804808 twice, TIC 240488136 once, TIC 2007804809 once).

## 2. Difference images (in transit vs out of transit)

`diff.py`: each event is measured as the in-transit mean (|dt| < 0.4 x 1.5 h) against a linear fit to the local out-of-transit window (0.75 x dur < |dt| < 0.3 d). Events are then stacked. Per-pixel noise comes from the robust spread of 120 stacks at random phases, with the stamp-wide common mode removed. The image is fitted with the TESS PRF from the `tess_prf` package that TESS-Localize uses. Figures: `out/diff_s15.png`, `out/diff_s16.png` (panels: out-of-transit image, difference image, difference/noise; green + marks the target, magenta x the best free-position fit).

| Sector | Events | Free-position PRF best fit | 1σ box (px) | Target-only PRF amplitude (needed for 366 ppm) |
|---|---|---|---|---|
| 15 | 12 | 13.3" from target | x 6.9–7.5, y 7.8–8.4 (target at 7.83, 7.96: inside y, 0.3 px outside x) | 6.89 ± 1.77 e/s (6.34 needed) |
| 16 | 12 | 42.8" from target, near G18.5 star …1772032 | x 4.4–5.8, y 6.2–7.3 (target at 7.25, 7.34: outside) | 3.35 ± 1.49 e/s (5.83 needed) |

- [verified] S15: the dimming sits on the target's pixels. Its amplitude matches the 366 ppm needed on the target (6.89 ± 1.77 against 6.34).
- [verified] S16: the target-only amplitude is 3.35 ± 1.49, which is 1.7σ below what is needed. The free fit's minimum is a single-pixel +2.5σ excess about 2 px west of the target, and no bright star sits there. Δχ² free-best vs target = 12.7. [inferred] With 12 events at about 2σ per pixel, a 2–3σ noise pixel can pull the free minimum; one sector at this level does not show an offset.
- [not done in budget] The 2-minute sectors (56, 76, 82, 83) were still running when this was written; they were slow because machine load was about 250–330. When they finish they will write to `out/diff_s{56,76,82,83}.{json,png,npz}`, and `post.py 56 76 82 83` will extend the neighbour table. Nothing below depends on them.
- Earlier SPOC flux-weighted centroid check (`hunt-ffi/centroid_240488148.txt`, verified): in-transit shifts of 0.1–6.3 mas, all |z| < 2.4. That is consistent with an on-target source.

## 3. Depth each neighbour would need, and whether the pixels allow it (S15 + S16)

Required depth = (aperture flux deficit) / (that star's PRF flux inside the SPOC aperture). The aperture deficit is 366 ppm x CROWDSAP x aperture flux, using the SPOC apertures from `hunt-ffi/lc/240488148/spoc_s*.fits`. The required amplitude A_req is the e/s drop that star would have to show at its own position. Δχ² is the forced fit at A_req (PRF at the star's Gaia position plus a constant) minus the target-only free fit, summed over S15 and S16. A star is ruled out at Δχ² > 25 (about 5σ). TRICERATOPS's own Gaussian-PSF required depths (`out/trc_depths.csv`) agree to within a few percent of the depth (for example, TIC 240488141 needs 5.1% by TRICERATOPS and 6.2% by PRF).

| TIC | Sep (") | Tmag | Depth needed | A_req S15 / S16 (e/s) | A measured S15 / S16 (e/s) | Δχ² | Verdict |
|---|---|---|---|---|---|---|---|
| 240488148 (target) | 0.0 | 10.17 | 0.037% | 6.3 / 5.8 | 6.89±1.77 / 3.35±1.49 | 5 | consistent |
| 2007804809 | 3.7 | 18.13 | 65% | 6.4 / 5.8 | 6.75±1.59 / 3.67±1.38 | −1 | **not excluded** |
| 2007804808 | 7.3 | 18.13 | 66% | 6.5 / 5.8 | 5.67±1.35 / 3.58±1.25 | 0 | **not excluded** |
| 240488136 | 11.2 | 18.22 | 74% | 6.5 / 6.1 | 4.08±1.50 / 2.20±1.42 | 23 | **not excluded** (marginal) |
| 240488160 | 22.5 | 17.30 | 34% | 7.5 / 6.1 | 2.71±0.70 / 1.92±0.71 | 80 | ruled out |
| 240488127 | 30.4 | 13.68 | 1.7% | 10.0 / 8.9 | 1.30±0.62 / 0.11±0.52 | 496 | ruled out |
| 2007804814 | 30.5 | 17.65 | 65% | 11.8 / 7.2 | 1.26±0.46 / 1.06±0.36 | 822 | ruled out |
| 240488129 | 30.8 | 15.31 | 7.2% | 10.0 / 8.0 | 1.42±0.61 / 0.12±0.48 | 496 | ruled out |
| 240488158 | 30.9 | 16.89 | 32% | 9.9 / 8.5 | 1.10±0.43 / −0.55±0.50 | 776 | ruled out |
| 240488124 | 31.4 | 17.42 | 53% | 9.2 / 9.9 | 0.36±0.72 / 0.66±0.72 | 336 | ruled out |
| 2007804784 | 31.5 | 18.03 | 90% | 10.3 / 8.0 | 1.37±0.58 / 0.10±0.46 | 554 | ruled out |
| 240488137 | 35.1 | 17.54 | 68% | 12.4 / 9.5 | 0.20±0.39 / −0.23±0.54 | 1333 | ruled out |
| 240488118 | 36.0 | 15.82 | 18% | 14.1 / 13.7 | 0.75±0.50 / 0.04±0.45 | 1652 | ruled out |
| 2007804793 | 36.8 | 17.73 | 95% | 11.3 / 14.5 | 0.31±0.58 / 0.65±0.52 | 1091 | ruled out |
| 2007804782 | 40.1 | 18.09 | >100% | 16.3 / 9.5 | 0.78±0.40 / −0.13±0.34 | 2349 | ruled out (also impossible) |
| **240488141** | **46.3** | **14.55** | **6.2%** | **20.9 / 10.2** | **1.02±0.34 / −0.05±0.24** | **5346** | **ruled out** |
| 240488110 | 47.9 | 13.51 | 3.7% | 23.1 / 25.9 | −0.20±0.44 / 0.68±0.38 | 7286 | ruled out |
| 240488120 | 48.2 | 16.99 | 88% | 29.6 / 17.4 | −0.08±0.40 / 0.51±0.38 | 7430 | ruled out |
| 240488163 | 48.7 | 17.18 | 61% | 17.8 / 9.7 | 0.58±0.38 / 0.63±0.26 | 3381 | ruled out |
| 240488130 | 51.8 | 14.74 | 16% | 47.0 / 22.0 | 0.05±0.29 / 0.33±0.42 | 29554 | ruled out |
| 240488128 | 52.5 | 17.35 | >100% | 41.3 / 16.2 | 0.58±0.28 / −0.23±0.20 | 29422 | ruled out |
| 240488167 | 53.3 | 17.24 | 87% | 25.3 / 11.4 | 0.37±0.36 / 0.81±0.28 | 6476 | ruled out |
| 240488151 | 54.4 | 14.58 | 18% | 53.8 / 34.9 | −0.12±0.32 / −0.28±0.35 | 39131 | ruled out |
| 240488099 | 56.3 | 16.35 | >100% | 63.1 / 64.4 | 0.15±0.30 / 0.29±0.22 | 136304 | ruled out |
| 240488181 | 56.3 | 15.80 | 97% | 101.0 / 53.4 | 0.44±0.16 / −0.13±0.18 | 492894 | ruled out |
| 240488119 | 57.4 | 14.50 | 20% | 74.7 / 31.4 | −0.11±0.35 / 0.42±0.33 | 56660 | ruled out |
| 240488134 | 58.1 | 13.26 | 7.7% | 87.9 / 40.2 | −0.19±0.33 / 0.11±0.44 | 82928 | ruled out |
| 240488095 | 62.5 | 14.91 | 70% | 157.2 / 96.2 | 0.27±0.28 / 0.12±0.22 | 513429 | ruled out |
| 240488107 | 64.7 | 15.69 | >100% | 145.4 / 69.1 | −0.09±0.34 / −0.25±0.24 | 274543 | ruled out |
| 240488149 | 68.9 | 13.11 | 17% | 235.7 / 87.9 | −0.17±0.39 / −0.28±0.36 | 448001 | ruled out |

These are the 30 stars where TRICERATOPS has tdepth > 0, meaning they get flux in the aperture. TRICERATOPS assigns no NEB/NTP probability to its other 768 stars, so they have no scenario to drop.

**Ruled out (26 stars, all ≥ 22.5"):** 240488160, 240488127, 2007804814, 240488129, 240488158, 240488124, 2007804784, 240488137, 240488118, 2007804793, 2007804782, 240488141, 240488110, 240488120, 240488163, 240488130, 240488128, 240488167, 240488151, 240488099, 240488181, 240488119, 240488134, 240488095, 240488107, 240488149 (`out/drop_ge20.txt`).

## 4. TRICERATOPS rerun with ruled-out neighbours removed

`trc_rerun.py`: same light curve (`hunt-ffi/lc/240488148_TESS-SPOC_masked.npz`), same binning, depth 366.1 ppm (dur 0.0625 d), same SPOC apertures and the same pickled target as `hunt-ffi/trc`. Stars are dropped with `target.remove_star()`, and `pix_coords` is filtered to match. Settings: N = 1e6, nsamples 20.

| Run | Stars dropped | FPP | NFPP | NFPP < 0.001? |
|---|---|---|---|---|
| original (hunt-ffi, 5 seeds) | 0 | 0.077–0.087 | 0.0035–0.0038 | no |
| baseline rerun here, seed 0 | 0 | 0.0868 | 0.00376 | no (reproduces original) |
| drop TIC 240488141 only, seed 0 | 1 | 0.0835 | 1.2e-4 | **yes** |
| drop all 26 ruled out, seed 0 | 26 | 0.0833 | 7.0e-9 | **yes** |
| seed 1 | 26 | 0.0760 | 7.5e-9 | **yes** |
| seed 2 | 26 | 0.0773 | 6.9e-9 | **yes** |

- [verified] After the drop, all of the remaining false-positive probability is on the target itself: PTP 0.143, STP 0.082, DTP 0.021, BTP 0.0005, BEB 0.0002, BEBx2P 0.0002 (`out/trc_probs_drop20_0.csv`). These are the unresolved-companion and background scenarios that pixel work cannot touch. They are what keeps FPP at about 0.08.
- [inferred] The NFPP bar of < 0.001 is met. The FPP bar of < 0.015 used for statistical validation is not met on TESS data alone. The FPP is dominated by PTP and STP, the "transit on an unresolved bound companion" scenarios. A contrast curve (`contrast_curve_file`) is the input that constrains those.

## Caveats

- [verified] The neighbour exclusions rest on the two 30-minute FFI sectors (24 events). The 2-minute-sector difference images had not finished. TESS-Localize did use all six sectors.
- [inferred] The forced-fit test uses `tess_prf` models at Gaia positions. A sub-pixel WCS error of about 0.1–0.2 px would not change any verdict for stars ≥ 22": their Δχ² margins are ≥ 3x the threshold.
- [guess] TIC 240488136 (11.2", Δχ² = 23) might cross 25 once the 2-minute sectors are added, but it would need a 74% eclipse on a G = 19 star. TRICERATOPS already gives it about 0.

## Files

- `work/tic240488148-localize/loc.py`, `out/loc_s*.json|csv|png` for TESS-Localize
- `work/tic240488148-localize/diff.py`, `out/diff_s15.*`, `out/diff_s16.*` for difference images
- `work/tic240488148-localize/post.py`, `out/neighbour_tests.csv`, `trc_depths.py`, `out/trc_depths.csv` for required depths and exclusion
- `work/tic240488148-localize/trc_rerun.py`, `out/trc_{base,drop141,drop20}_*.json`, `out/trc_probs_*.csv` for TRICERATOPS

## Addendum: Sector 83 (2-minute FFI, 15 events), finished after the main write-up

`out/diff_s83.{json,png,npz}`, `out/neighbour_tests_s83.csv` [verified]
- Target-only PRF amplitude 3.91 ± 1.12 e/s, against 5.77 needed for 366 ppm (1.7σ low, same as S16). 3x3 aperture depth 269 ± 58 ppm. The free-position best fit is 13.6" from the target (x 7.9–8.3, y 6.6–7.1; target at 7.96, 7.57, just outside the 1σ box in y).
- Every one of the 26 stars ruled out in S15+S16 is independently ruled out by S83 alone (Δχ² 98 to 27,846). TIC 240488141 needs 8.94 e/s, measures −0.02 ± 0.35, Δχ² = 699.
- The three close faint stars give S83 Δχ² of 6.5 (TIC 2007804809), 10.8 (TIC 2007804808) and 17.0 (TIC 240488136). Summed with S15+S16 that is 5.9, 10.3 and **40.0**. [inferred] TIC 240488136 (11.2") now passes the Δχ² > 25 cut over three sectors. The other two remain not excluded.
- Sectors 56, 76 and 82 were still running at this point (`out/diff_5*.log` etc.).

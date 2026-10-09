# TIC 240488148.01: a 1.4 Earth-radius, 1.51-day transiting planet candidate

Discovered by Advait Paliwal on 2026-10-08 in public TESS full-frame-image photometry.

## Summary

TIC 240488148 (TYC 3191-503-1) is a K dwarf (Teff 4403 K, 0.67 R☉) at 58.7 pc. TESS-SPOC full-frame-image light curves from Sectors 15, 16, 56, 76, 82 and 83 (2019 to 2024) show a transit-like signal with period 1.5125149 ± 0.0000029 d, depth 411 ± 28 ppm and duration 1.41 h, which implies a radius of 1.37 ± 0.14 Earth radii. The signal is present in every sector and in two independent light-curve pipelines (TESS-SPOC and QLP), odd and even depths agree, and there is no secondary eclipse. Difference images and TESS-Localize exclude every catalogued neighbour beyond 11″ as the source. As of 2026-10-09 the star is not a TOI or a CTOI and is not in the NASA Exoplanet Archive. The TRICERATOPS false-positive probability is 0.081 without companion constraints and 0.024 with Gaia-based companion limits from MOLUSC, above the 0.015 threshold for statistical validation. **This is a planet candidate, not a validated or confirmed planet.** High-resolution imaging is needed to rule out unresolved companions.

| | |
|---|---|
| Star | TIC 240488148, RA 323.62829°, Dec +42.64403° (J2000), Tmag 10.17 |
| Period | 1.5125149 ± 0.0000029 d |
| Mid-transit (T0) | BJD_TDB 2460393.34348 ± 0.0013 |
| Depth | 411 ± 28 ppm |
| Duration | 1.41 h |
| Radius | 1.37 (+0.14 / −0.13) R⊕ |
| Equilibrium temperature | about 1130 K |
| False-positive probability | 0.024 (Gaia companion limits), 0.081 (none) |

## Contents

- `note/TIC240488148_note.md`, `note/TIC240488148_note.pdf`: the research note (abstract, star, data, signal, neighbour tests, false-positive analysis).
- `note/figure1_TIC240488148.png`: all six sectors folded on the fitted ephemeris, with the transit model and residuals.
- `ephemeris.csv`: predicted mid-transit times from 2026-10-09 to 2027-12-31.
- `code/`: the scripts used (search, vetting, transit fit, pixel localization, companion limits and false-positive runs, ephemeris).
- `data/`: small outputs of those scripts (fit posterior and results, vetting reports, difference images, TRICERATOPS and MOLUSC outputs).
- `data/exofop_ctoi_row.txt`: the prepared ExoFOP community-candidate parameter row (not yet uploaded).

## Next predicted transits

From `ephemeris.csv` (1σ timing uncertainty about 4 minutes now, growing to 5 minutes by end of 2027):

| Epoch | Mid-transit BJD_TDB | Mid-transit UTC | ±1σ |
|---|---|---|---|
| 615 | 2461323.54014 | 2026-10-10 00:52 | 4.0 min |
| 616 | 2461325.05266 | 2026-10-11 13:10 | 4.0 min |
| 617 | 2461326.56517 | 2026-10-13 01:28 | 4.0 min |
| 618 | 2461328.07769 | 2026-10-14 13:46 | 4.0 min |
| 619 | 2461329.59020 | 2026-10-16 02:04 | 4.0 min |
| 620 | 2461331.10272 | 2026-10-17 14:22 | 4.0 min |

Times are propagated from every posterior sample of the transit fit, so the correlation between T0 and P is kept. UTC times are geocentric (BJD_TDB minus the light-travel delay to the star, then TDB to UTC); ingress and egress are mid-transit ± 0.705 h. Epochs count from T0 = BJD_TDB 2460393.34348.

## Reproducing

Python 3.12 with numpy, scipy, astropy, lightkurve, transitleastsquares, wotan, batman-package, emcee, corner, tess-localize and triceratops. Package versions for the search and false-positive environments are in `data/companions/hunt-ffi-freeze.txt` and `data/companions/trc-freeze.txt`.

1. Light curves. Download the TESS-SPOC FFI light curves for TIC 240488148, Sectors 15, 16, 56, 76, 82 and 83, from MAST (`code/search/run1.py` does this and runs Transit Least Squares). FITS light curves and TESScut cutouts are not included here; they are public at MAST.
2. Transit fit. `TIC_LC_DIR=<folder with spoc_sNN.fits> python code/fit/fit.py` from inside `data/fit/` (it reads `_sources/` for the stellar and limb-darkening priors), then `python code/fit/post.py` for the derived parameters, corner plot and folded plot. The posterior used here is `data/fit/chain.npz`; the results are in `data/fit/result.json` and `data/fit/RESULT.md`.
3. Ephemeris. `python code/ephemeris.py data/fit/chain.npz ephemeris.csv`.
4. Vetting and neighbours. `code/search/` (vet, centroid, shape, xcheck, Gaia flags) and `code/localize/` (difference images, TESS-Localize, TRICERATOPS with excluded neighbours). Outputs are in `data/search/` and `data/localize/`.
5. Companion limits. `code/companions/` runs MOLUSC (https://github.com/woodml/MOLUSC, not bundled) on Gaia constraints and feeds the result to TRICERATOPS. Outputs are in `data/companions/`. The 1e6-row MOLUSC file passed to TRICERATOPS (`molusc_for_trc_N1e6.csv`) is not included because of its size; `code/companions/molusc_summary.py` regenerates it from the MOLUSC output.

Some scripts use paths relative to the original working folders (for example `work/hunt-ffi/`); adjust them to your layout.

## Credits

This work uses data collected by the TESS mission, funded by the NASA Science Mission Directorate, obtained from the Mikulski Archive for Space Telescopes (MAST), including TESS-SPOC (Caldwell et al. 2020) and QLP (Huang et al. 2020) light curves and TESScut cutouts. It uses data from the European Space Agency mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium, and the TESS Input Catalog (Stassun et al. 2019). Full references are in the note.

## Licence and citation

Code is under the MIT licence (`LICENSE`). Text, figures and data are under CC BY 4.0 (`LICENSE-CC-BY-4.0`). Please cite with `CITATION.cff`.

## Contact

Advait Paliwal, advaitpaliwal@gmail.com

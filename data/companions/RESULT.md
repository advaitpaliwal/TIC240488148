# TIC 240488148 (TYC 3191-503-1): bound-companion limits from archival data only, and TRICERATOPS re-run

Date: 2026-10-08. Labels: **verified** = read from a saved file or command output in this folder; **inferred** = my calculation or reading; **guess** = not checked.

## Bottom line

- **No archival high-resolution imaging, spectroscopy or ground RVs exist** for this star in any archive checked. The only spectroscopic data are the Gaia DR2/DR3 RVS radial velocities. (verified)
- **MOLUSC (Gaia RUWE 1.039 + Gaia imaging):** 41% of simulated stellar companions survive. With a 41-46% companion prior, the chance that an undetected bound stellar companion exists is **about 22-26%**, falling to about 14-17% if the constant Gaia RVS velocity is also used (my add-on test, not part of MOLUSC). (verified run output; the prior and the RV test are inferred)
- **TRICERATOPS with the Gaia-only constraint does not reach the validation bar.** FPP drops from 0.081 to about 0.024 (MOLUSC prior) and 0.062 (Gaia contrast curve), both above 0.015. **NFPP is 0.0037-0.0040, also above the 0.001 bar,** and companion limits cannot change it. (verified, numbers below)
- **The literature does not accept Gaia-only limits in place of imaging for a planet this small.** The 2026 TRICERATOPS imaging study found 22 of 22 validated planets below 1.7 Earth radii needed their imaging contrast curve to pass. (verified quote below)

## (a) Archival data on the star

Source: Gaia DR3 1967568881023980032, G 10.88, parallax 16.996 ± 0.013 mas (58.8 pc). (verified, `_sources/vizier_gaiadr3_0.csv`)

| Archive | Result | Evidence |
|---|---|---|
| ExoFOP | 0 imaging, 0 spectroscopy, 0 time series, 0 stellar companions, 0 files, no TOI/CTOI | `_sources/exofop_target_240488148.json` (verified) |
| SIMBAD | 3 papers, none observational: Gaia red-clump calibration (2018), Gaia Catalogue of Nearby Stars (2021), Bioverse ELT simulation (2023). No spectral type. | `_sources/simbad_tyc3191-503-1.txt`, `_sources/simbad_refs.csv` (verified) |
| Gaia DR3 non-single-star tables (two-body orbits, SB1, SB2, acceleration 7/9-param, VIM) | No entry; `non_single_star` = 0 | `q_archive.py` output, `vizier_gaiadr3_0.csv` (verified) |
| LAMOST DR5 LRS, DR7, DR8 MRS | No match within 5" | `q_archive.py` output (verified) |
| RAVE DR6, APOGEE DR17, GALAH DR3 | No match | same (verified) |
| WDS, El-Badry 2021 wide binaries, Tokovinin MSC | No match | same (verified) |
| Gaia RVS | DR3 RV −44.12 ± 0.18 km/s from 32 transits over 22 visibility periods; constancy p-value 0.83; robust amplitude 2.58 km/s. DR2 RV −43.90 ± 0.36. | `vizier_gaiadr3_0.csv`, SIMBAD (verified) |
| Gaia close-pair flags | RUWE 1.039, ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.022, duplicated_source 0 | `vizier_gaiadr3_0.csv` (verified) |
| Gaia astrometric excess noise | 0.088 mas, significance 12.1 | same (verified). Bright stars near G 11 commonly show a small significant excess noise from calibration; RUWE is the standard binarity test (inferred) |
| Gaia neighbours within 30" | 13 sources, all G ≥ 14.5 and parallax ≤ 2.1 mas, none co-moving | `vizier_gaiadr3_0.csv` (verified) |

Not checked: Keck KOA, ESO, Gemini and CFHT raw archives directly (guess: nothing, since ExoFOP and SIMBAD list nothing).

## (b) MOLUSC

Run: MOLUSC (github.com/woodml/MOLUSC, commit b26efbf) command-line mode, `--ruwe --gaia`, mass 0.69 Msun, default age 5 Gyr, 4 × 50,000 companions with seeds 11-14. Two changes to make it run, both in `molusc_run.py`: Gaia queries go to VizieR I/355/gaiadr3 (the ESA archive hangs), and a command-line bug (`star_jitter` never set) is patched. (verified)

| Quantity | Value | Label |
|---|---|---|
| Rejected by RUWE test | 18.3% of all simulated | verified |
| Rejected by Gaia imaging test | 34.7% | verified |
| Survive (all companions incl. substellar) | 48.1% (seeds: 48.05, 48.14, 48.12, 48.09%) | verified |
| Survive (stellar, m2 ≥ 0.08 Msun) | 41.3% | verified |
| Survivors with a(1−e) > 10 AU, the subset TRICERATOPS uses | 15.0% of all simulated | verified |
| P(undetected stellar companion), prior 41-46% | 22-26% | inferred (prior is a literature figure for K dwarfs, guess-level, not saved) |
| Plus Gaia RVS constancy (K > 3 × 0.64 km/s and P < 1038 d rejected) | 23.3% of stellar survive; P = 14-17% | inferred (`rv_filter.py`; per-epoch scatter back-computed from e_RV) |

Where companions survive (fraction rejected, `molusc_map.py`, verified):

| Projected separation | ΔG < 3 | ΔG 3-5 | ΔG 5-7 |
|---|---|---|---|
| < 1.1 AU (< 19 mas) | 0% | 0% | 0% |
| 1.1-5 AU | 24-48% | 2-10% | 1% |
| 5-20 AU | 98-100% | 39-81% | 8% |
| 20-60 AU | 87-99% | 37-66% | 12% |
| 60-130 AU (1-2.2") | 72-84% | 44-62% | 20% |

So Gaia removes most bright companions from 5 to 130 AU, but faint ones (ΔG > 5) and anything inside 5 AU are barely touched. (inferred from the table)

## (c) Gaia-derived contrast-vs-separation limit

`cc_gaia_dr2completeness_50pct_G.csv` and `_90pct_G.csv` are MOLUSC's bundled Gaia DR2 recovery curves (Ziegler et al. 2018 inside 4", Brandeker & Cataldi 2019 outside), the same file MOLUSC uses. (verified, `molusc/code/gaia_contrast.txt`; source description verified in `_sources/wood2021.txt` lines 510-545)

50% recovery, ΔG: 0" at ≤ 0.7", 1.9 at 1.0", 4.9 at 2.0", 6.3 at 3.0", 8.9 at 4.2", 13.0 at 10". The 90% curve starts at 1.7".

Label: this is a **statistical recovery rate for Gaia DR2 binaries in general, not a measured detection limit around this star**. It is not equivalent to a speckle or AO curve: speckle reaches ΔG ≈ 4-5 at 0.1-0.2", where this curve gives 0. Against TRICERATOPS's default (everything unresolved inside 2.2", nothing beyond), it is tighter for ΔG < 5 and looser for ΔG > 5. (inferred)

## (d) TRICERATOPS re-run

Setup: TRICERATOPS 1.1.0, same light curve, binning, aperture pickle and Gaia field-star file as `../hunt-ffi/trc/`, N = 1e6, 30-min exposure. 10 runs per mode, seeds 1000-1009. One upstream bug patched in `trc/run_fpp.py`: with pandas 3, the `molusc_file` branch fails on a read-only array. (verified)

The light-curve file `../hunt-ffi/lc/240488148_TESS-SPOC_masked.npz` gives a depth of 394.7 ppm with the default `dur` = 0.06 d, versus 366.1 ppm in the earlier run, with an identical out-of-transit scatter (70.42 ppm). `dur` = 0.065-0.07 d reproduces about 360-366 ppm, so the earlier run probably passed a different `dur` (inferred). So I re-ran a no-constraint baseline on today's depth. (verified, `res_*.json`)

| Mode | Runs | FPP mean ± SD | NFPP mean ± SD | TP | PTP | STP | DTP | vs FPP < 0.015 |
|---|---|---|---|---|---|---|---|---|
| No constraint (TRICERATOPS default 2.2") | 5 | 0.0812 ± 0.0034 | 0.00365 ± 0.00011 | 0.760 | 0.138 | 0.077 | 0.021 | fails (5.4×) |
| Gaia 50%-recovery contrast curve (`filt='Vis'`) | 10 | 0.0620 ± 0.0020 | 0.00367 ± 0.00012 | 0.750 | 0.113 | 0.058 | 0.075 | fails (4.1×) |
| MOLUSC posterior (`molusc_file`, RUWE + Gaia imaging) | 10 | 0.0235 ± 0.0035 | 0.00403 ± 0.00014 | 0.825 | 0.128 | 0.019 | 0.023 | fails (1.6×); best run 0.0198 |

Earlier hunt-ffi runs at 366 ppm depth: FPP 0.081 ± 0.003 over 5 runs (`../hunt-ffi/trc/out/res_240488148_*.json`), matching the baseline here. Per-run values are in `trc/out/res_*.json`; run `summarize.py` to reproduce. (verified)

Bar: FPP < 0.015 and NFPP < 0.001 (Giacalone et al. 2021, verified `_sources/giacalone2021.txt` line 91).

What the numbers show (inferred from the table):
- **MOLUSC mode** removes most of the companion-hosted scenario (STP 0.077 → 0.019). The diluted-planet scenario PTP stays at 0.128 because TRICERATOPS's PTP still allows faint companions with ΔG > 5 inside a few AU, and Gaia cannot reach those. All 10 runs fall between 0.0198 and 0.0289; none is under 0.015.
- **The Gaia contrast curve** helps less than MOLUSC and moves probability into DTP (0.021 → 0.075). The curve is weaker than the 2.2" default for faint stars at 2-4", so more unresolved background stars are allowed.
- **NFPP 0.0037-0.0040 fails the < 0.001 bar in every mode.** It comes from nearby star TIC 240488141 (NEB + NTP, about 0.004; `trc/out/probs_*.csv`). Companion limits cannot fix this; only on-target ground photometry or a pixel-level localization can.
- So even with the best Gaia-only input, this candidate is a "likely planet" (FPP < 0.5) but not a "validated planet" by either criterion.

## (e) Do validation papers accept Gaia-only limits instead of imaging?

**No paper I found validates a sub-2-Earth-radius planet with TRICERATOPS on Gaia-only companion limits.** The evidence points the other way:

1. Collier, Giacalone et al. 2026 (arXiv 2608.21560, accepted AAS Journals; `_sources/collier2026.txt` lines 2779-2785): "Without contrast curve constraints, 72% of our validated planets with adopted contrast curves would fail the FPP<0.015 threshold ... 100% of validated planets below 1.7R⊕ require the contrast curve to achieve validation". Lines 2489-2492: "Below 1.7R ⊕, 22 of 22 VPs with contrast curves (100%) are CC-dependent". (verified) They use RUWE ≤ 1.4 only as an extra quality cut, not as a replacement (lines 934-938). (verified)
2. Same paper, lines 2822-2825: the needed property is "sensitivity to faint companions close to the star, with Gaia astrometry independently flagging a fraction of close unresolved companions". "A fraction" means Gaia only covers part of that space. (verified)
3. Giacalone et al. 2021 (`_sources/giacalone2021.txt` lines 1056-1067): the companion weight uses "the contrast curve obtained from high-resolution imaging"; with no imaging "the angular separation ... is assumed to be 2.′′2 (Brown et al. 2018)". So the method's own Gaia-based default is already built in. (verified)
4. Wood, Mann & Kraus 2021 (MOLUSC; `_sources/wood2021.txt` lines 1152-1156): "If only constraints from Gaia are sufficient for false-positive analysis, they may make additional HRI observations unnecessary". The authors frame this as conditional. They also write (lines 1143-1147): "Even with full follow-up data (RV and HRI measurements), we are often unable to completely rule out a stellar-mass companion". (verified)
5. The TRICERATOPS author recommends MOLUSC output as an input: "While it is optional, I recommend people use it" (`_sources/triceratops_TSCIII_tutorial.ipynb`). His example (TOI-1228) pairs it with a speckle contrast curve, not in place of one. (verified)
6. VaTEST II (Mistry et al. 2023, AJ, from search result text, not saved as a file): RUWE ≥ 1.4 is used only to remove targets, and validation used "high-resolution imaging, ground-based photometry, and transit photometry data". (verified from search highlight only)

Verdict (inferred): a Gaia/MOLUSC limit is a legitimate *input* that published work uses alongside imaging. Using it *instead of* imaging to claim validation of a 1.4-Earth-radius planet would sit outside current practice, and here it does not clear the bar anyway.

## Files

- `q_archive.py`, `_sources/vizier_*.csv`: archive sweep
- `molusc_run.py`, `molusc_out/m50k_s1{1..4}_*.csv`, `molusc_summary.py`, `molusc_map.py`, `rv_filter.py`: MOLUSC
- `molusc_out/molusc_for_trc_N1e6.csv`: MOLUSC survivors resampled to 48.1% of 1e6 rows. TRICERATOPS pads the rest with no companion; this follows the upstream example file. It treats every star as having a companion a priori, which is conservative (inferred)
- `cc_gaia_dr2completeness_{50,90}pct_G.csv`: Gaia contrast curves
- `trc/run_fpp.py`, `trc/out/res_*.json`, `trc/out/probs_*.csv`, `summarize.py`: TRICERATOPS
- `trc/out/discarded_clipped_cc/`: one early run with a wrongly built curve (clipped at G = 18 with a 60" tail, which inflated DTP to 0.40); excluded
- `_sources/`: ExoFOP, SIMBAD, Giacalone 2021, Collier 2026, Wood 2021 PDFs and text with `.url` notes

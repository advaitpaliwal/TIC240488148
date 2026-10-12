# A 1.4 Earth-radius Planet Candidate on a 1.51-day Orbit around the K Dwarf TIC 240488148

Advait Paliwal
Independent researcher, Sterling Heights, MI 48314, USA
advaitpaliwal@gmail.com

## Abstract

I report a transit-like signal on TIC 240488148 (TYC 3191-503-1), a K dwarf at 58.7 pc, in TESS full-frame-image photometry from six sectors spanning 2019–2024. A transit fit gives a period of 1.5125149 ± 0.0000029 d, a depth of 411 ± 28 ppm and a duration of 1.41 h, implying a radius of 1.37 ± 0.14 R⊕. The star is not a TOI or CTOI, and no SPOC report exists for it. The signal appears in every sector and in two independent pipelines, with consistent odd and even depths and no secondary eclipse. Difference images in all six sectors rule out the brightest nearby star and every catalogued neighbour beyond 20″ as the source. With Gaia companion limits and no neighbours removed, TRICERATOPS gives FPP = 0.024 ± 0.004, above the 0.015 validation threshold, so the candidate awaits high-resolution imaging.

## Star and search

TIC 240488148 (Stassun et al. 2019) has T = 10.17, Teff = 4403 ± 127 K, R★ = 0.672 ± 0.057 R☉ and M★ = 0.690 ± 0.083 M☉. Its Gaia DR3 RUWE is 1.039, and its Gaia radial velocity, −44.12 ± 0.18 km s⁻¹, is constant over 32 epochs (p = 0.83; Gaia Collaboration et al. 2023). ExoFOP lists no imaging, spectroscopy or radial-velocity follow-up, and the star is absent from the Gaia non-single-star tables and the wide-binary catalogue of El-Badry et al. (2021).

The star was never observed at 2-minute cadence before Sector 56 and has no SPOC data-validation report. I searched the TESS-SPOC full-frame-image light curves (Caldwell et al. 2020) of 1,096 nearby K and M dwarfs lacking 2-minute data with Transit Least Squares (Hippke & Heller 2019) over 0.5–12 d. Before the search, each light curve was flattened with the lightkurve Savitzky–Golay filter (window 0.75 d, no transit mask) and binned to 30 minutes. The same pipeline recovers TOI-270 at its known period. TIC 240488148 gave SDE 29.2 from Sectors 15 and 16 (30-minute cadence) and 56, 76, 82 and 83 (200-second cadence).

## Transit fit

For the fit I detrended each sector with wotan (Hippke et al. 2019) using a biweight filter with a 0.75 d window. Points within ±1.875 h of each predicted transit were removed before detrending, and the trend was interpolated across them. I then fitted 10,263 points within 4.5 h of transit with batman (Kreidberg 2015) and emcee (Foreman-Mackey et al. 2013). The model is supersampled to each cadence, assumes a circular orbit, and uses a stellar-density prior from the TIC and quadratic limb darkening from Claret (2017). The fit gives P = 1.5125149 ± 0.0000029 d, T0 = 2460393.3435 (+0.0015/−0.0012) BJD_TDB, Rp/R★ = 0.0187 ± 0.0009, a/R★ = 7.6 (+0.5/−0.8), b = 0.42 (+0.18/−0.24), a depth of 411 ± 28 ppm and a duration of 1.41 (+0.07/−0.06) h (Figure 1). This implies Rp = 1.37 (+0.14/−0.13) R⊕, a = 0.0228 ± 0.0009 AU, an insolation of 294 (+32/−27) S⊕ and an equilibrium temperature of 1133 (+68/−51) K. The stellar density implied by the transit, 3.6 ± 0.9 g cm⁻³, agrees with the TIC value of 3.2 ± 1.0 g cm⁻³.

## Vetting

Box-fit depths per sector are 337 ± 91, 230 ± 73, 413 ± 57, 438 ± 90, 433 ± 59 and 340 ± 81 ppm. Each is the median out-of-transit flux minus the mean in-transit flux, and each error is the scatter of the same statistic at 300 random phases. Against their weighted mean (374 ± 29 ppm) they give χ² = 6.2 for 5 degrees of freedom. The independent QLP light curves (Huang et al. 2020a, 2020b; Kunimoto et al. 2021) recover the signal in all five sectors they cover, at 412 ± 33 ppm (χ² = 9.6 for 4 degrees of freedom). Odd and even transits agree to 0.67σ, the depth at phase 0.5 is 30 ± 23 ppm, and folding at half and twice the period confirms the stated period. A single sector shows individual transits at S/N ≈ 2 (Figure 1, lower panels).

I made in-transit difference images from TESScut cutouts (Brasseur et al. 2019) for all six sectors. In each I fitted the TESS PRF at the target and at each catalogued neighbour, scaled to the flux loss that neighbour would need to produce the signal. The quoted Δχ² is the χ² of that forced neighbour model minus that of a free target model, over the 15 × 15 pixel stamp, with per-pixel noise from 120 difference images at random phases. The flux deficit sits on the target. TIC 240488141 (T = 14.6, 46″ away) would need a 5.5% eclipse; placing the signal on it is ruled out in every sector, by Δχ² = 699–21,069 per sector (42,124 in total). The 26 other catalogued stars between 22.5″ and 69″ are each ruled out in every sector. Three faint stars within 12″ (G ≈ 19–20) cannot be separated from the target at TESS resolution, but each would need an eclipse of about 70%.

## False-positive probability

I ran TRICERATOPS (Giacalone & Dressing 2020; Giacalone et al. 2021) ten times with 10⁶ draws each. Companion limits came from MOLUSC (Wood et al. 2021), using its Gaia RUWE and contrast tests. It generated 200,000 companions, kept 48%, and the kept companions replaced TRICERATOPS's default companion prior for the bound-companion scenarios. With all 798 Gaia sources retained, the result is FPP = 0.024 ± 0.004 and NFPP = 0.004. Of that NFPP, 97% comes from TIC 240488141. Following Giacalone et al. (2021), who omit nearby-star scenarios once ruled out, removing only this star gives NFPP = 1.3 × 10⁻⁴ and FPP = 0.020 ± 0.004. This exclusion rests on the TESS pixels, not on ground-based photometry. Either way the FPP is above 0.015 and is dominated by a planet orbiting an unresolved bound companion, so I make no validation claim (cf. Collier et al. 2026).

## Summary

TIC 240488148.01 is a 1.37 R⊕ candidate on an ultra-short-period orbit around a bright, nearby K dwarf. Speckle or adaptive-optics imaging can test the remaining companion scenarios, ground-based photometry on a 1 m telescope can confirm the transit is on target, and precise radial velocities could measure its mass (K ≈ 1.5–2 m s⁻¹ for a rocky composition). Code, data and predicted transit times are at https://github.com/advaitpaliwal/TIC240488148.

## Acknowledgments

This paper includes data collected by the TESS mission, funded by the NASA Science Mission Directorate, and obtained from MAST. It uses data from the ESA mission Gaia, processed by the Gaia DPAC.

## Figure 1

TESS photometry of TIC 240488148. Top: TESS-SPOC data from Sectors 15, 16, 56, 76, 82 and 83 folded on the fitted ephemeris (grey), with 15-minute bins (black) and the median model (orange). Middle: the unfolded Sector 82 light curve (200 s), PDCSAP flux with the wotan trend (offset) above the detrended flux in 30-minute bins, with predicted transit times marked. Bottom: four evenly spaced Sector 82 transits in 30-minute bins.

## References

Brasseur, C. E., Phillip, C., Fleming, S. W., Mullally, S. E., & White, R. L. 2019, Astrocut: Tools for creating cutouts of TESS images, ascl:1905.007
Caldwell, D. A., Tenenbaum, P., Twicken, J. D., et al. 2020, RNAAS, 4, 201
Claret, A. 2017, A&A, 600, A30
Collier, M., Giacalone, S., Derfer, B., et al. 2026, AAS Journals, accepted, arXiv:2608.21560
El-Badry, K., Rix, H.-W., & Heintz, T. M. 2021, MNRAS, 506, 2269
Foreman-Mackey, D., Hogg, D. W., Lang, D., & Goodman, J. 2013, PASP, 125, 306
Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. 2023, A&A, 674, A1
Giacalone, S., & Dressing, C. D. 2020, triceratops, ascl:2002.004
Giacalone, S., Dressing, C. D., Jensen, E. L. N., et al. 2021, AJ, 161, 24
Hippke, M., David, T. J., Mulders, G. D., & Heller, R. 2019, AJ, 158, 143
Hippke, M., & Heller, R. 2019, A&A, 623, A39
Huang, C. X., Vanderburg, A., Pál, A., et al. 2020a, RNAAS, 4, 204
Huang, C. X., Vanderburg, A., Pál, A., et al. 2020b, RNAAS, 4, 206
Kreidberg, L. 2015, PASP, 127, 1161
Kunimoto, M., Huang, C., Tey, E., et al. 2021, RNAAS, 5, 234
Stassun, K. G., Oelkers, R. J., Paegert, M., et al. 2019, AJ, 158, 138
Wood, M. L., Mann, A. W., & Kraus, A. L. 2021, AJ, 161, 85

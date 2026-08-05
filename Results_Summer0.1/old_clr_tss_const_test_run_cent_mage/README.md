# [DEPRECATED] Old CLR/TSS Constant Test Run (cent vs mage)

**Status: SUPERSEDED** - replaced by `coeff_sm_single_run_test_cent_mage/` and `coeff_sm_regression_analysis_all_pairs/` which use the correct `coefficient.sm` function from ALDEx3 v1.2.0.

## What was wrong

This was a test run of the old custom `clr_const` and `tss_const` functions before the full analysis. Same bug as the full runs: both functions applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope) instead of only the slope. The intercept should remain at 0 with no uncertainty.

The correct approach is `coefficient.sm` with `c.mu = c(0, k)` and `c.cor = diag(c(0, gamma^2))`.

## Script

`summary_focused_analysis_2.0_TEST.R` (copy included; source in `code/`)

## Method (old, incorrect)

- Custom `clr_const` and `tss_const` functions applied k and gamma to all P coefficients
- Single gamma: 0.3
- Single comparison: cent vs mage only
- Both CLR_constant/ and TSS_constant/ subdirectories
- Dataset: Bian et al. 2017 (1117 OTUs)

## Output structure

- `CLR_constant/gamma_0.3/{data,plots}/` - CLR side results
- `TSS_constant/gamma_0.3/{data,plots}/` - TSS side results

## Replacement

Use `coeff_sm_single_run_test_cent_mage/` for single-run testing or `coeff_sm_regression_analysis_all_pairs/` for full analysis.

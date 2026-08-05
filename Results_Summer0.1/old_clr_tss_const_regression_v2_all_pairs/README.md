# [DEPRECATED] Old CLR/TSS Constant Regression v2.0 (all pairs)

**Status: SUPERSEDED** - replaced by `coeff_sm_regression_analysis_all_pairs/` which uses the correct `coefficient.sm` function from ALDEx3 v1.2.0.

This is the latest version of the old custom constant analysis (v2.0). Earlier versions exist in `old_clr_const_regression_analysis_all_pairs/` and `old_tss_const_regression_analysis_all_pairs/` but are older runs.

## What was wrong

Used custom `clr_const` and `tss_const` functions that applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope) via `rnorm(P * ns, k, gamma)`. The intercept should remain at 0 with no uncertainty - only the slope should be perturbed. The correct approach is `coefficient.sm` with `c.mu = c(0, k)` and `c.cor = diag(c(0, gamma^2))`.

## Script

`summary_focused_analysis_2.0.R` (copy included; source in `code/`)

## Method (old, incorrect)

- Custom `clr_const` and `tss_const` functions applied k and gamma to all P coefficients
- Two-part analysis: CLR(const) vs TSS(std), then CLR(std) vs TSS(const)
- Gammas: 0.1, 0.3, 0.5
- Constants: -1 to 1 by 0.1
- Both comparisons: cent vs mage, eld vs kin
- Dataset: Bian et al. 2017 (1117 OTUs)

## Output structure

- `pval_analysis0.2CLR_constant/gamma_{0.1,0.3,0.5}/{data,plots}/`
- `pval_analysis0.2TSS_constant/gamma_{0.1,0.3,0.5}/{data,plots}/`

## Replacement

Use `coeff_sm_regression_analysis_all_pairs/` instead, which runs `coefficient.sm` (ALDEx3 v1.2.0) with proper three-way comparison (coeff vs CLR vs TSS) and the unified 5-category system.

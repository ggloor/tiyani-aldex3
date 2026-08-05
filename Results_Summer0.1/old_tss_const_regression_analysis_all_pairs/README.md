# [DEPRECATED] Old TSS Constant Regression Analysis (all pairs)

**Status: SUPERSEDED** - replaced by `coeff_sm_regression_analysis_all_pairs/` which uses the correct `coefficient.sm` function from ALDEx3 v1.2.0.

## What was wrong

This analysis used a custom `tss_const` function that manually applied a constant k and uncertainty gamma to the TSS scale model. The problem: `tss_const` applied the same mean k and SD gamma to ALL design matrix coefficients (both intercept and slope). This is incorrect because:

- The intercept should remain unperturbed (mean = 0, no uncertainty)
- Only the slope (treatment effect) should receive the perturbation

The correct approach is `coefficient.sm` with `c.mu = c(0, k)` and `c.cor = diag(c(0, gamma^2))`, which correctly targets only the slope coefficient. This was implemented in `summary_focused_analysis_3.0.R`.

## Script

`summary_focused_analysis_2.0.R` (copy included; source in `code/`)

## Method (old, incorrect)

- Custom `tss_const` function applied k and gamma to all P coefficients
- Compared TSS(constant) vs standard CLR
- Two-part analysis: TSS side and CLR side separately
- Gammas: 0.1, 0.3, 0.5
- Both comparisons: cent vs mage, eld vs kin
- Dataset: Bian et al. 2017 (1117 OTUs)

## Output structure

- `gamma_{0.1,0.3,0.5}/data/` - summary stats CSV, red dot details CSV, Rda files
- `gamma_{0.1,0.3,0.5}/plots/` - regression PDF, scatter PDF, effect PDF

## Replacement

Use `coeff_sm_regression_analysis_all_pairs/` instead, which runs `coefficient.sm` (ALDEx3 v1.2.0) with proper three-way comparison (coeff vs CLR vs TSS) and the unified 5-category system.

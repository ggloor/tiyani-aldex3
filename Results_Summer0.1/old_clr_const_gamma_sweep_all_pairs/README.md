# [DEPRECATED] Old CLR Constant Gamma Sweep (all 21 pairs)

**Status: SUPERSEDED** - replaced by `coeff_sm_effect_sweep_cent_mage/` which uses the correct `coefficient.sm` function from ALDEx3 v1.2.0.

## What was wrong

All scripts used old custom `clr_const` or `tss_const` functions that applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope) via `rnorm(P * ns, k, gamma)`. The intercept should remain at 0 - only the slope should be perturbed. Correct approach: `coefficient.sm` with `c.mu = c(0, k)` and `c.cor = diag(c(0, gamma^2))`.

## Folder structure

- `scripts/` - all R scripts that produced these results
- `data/` - Rda regression data files (cent vs mage, eld vs kin, CLR/TSS sides)
- `pval_scatter_plots/` - -log10(FDR) scatter plots across gamma and constant combinations
- `effect_cross_plots/` - effect cross plots (estimate vs std error) across constants

See subfolder READMEs for details.

## Replacement

Use `coeff_sm_effect_sweep_cent_mage/` for effect sweep analysis with the correct coefficient.sm function.

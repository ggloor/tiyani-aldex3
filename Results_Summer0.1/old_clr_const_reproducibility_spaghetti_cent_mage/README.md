# [DEPRECATED] Old CLR Constant Reproducibility Spaghetti Plots (cent vs mage)

**Status: SUPERSEDED** - replaced by `coeff_sm_reproducibility_spaghetti_cent_mage/` which uses the correct `coefficient.sm` function from ALDEx3 v1.2.0 and compares all three methods (coeff vs CLR vs TSS).

## What was wrong

1. **Incorrect custom function**: Used `clr_const` which applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope). Only the slope should be perturbed.
2. **Two-way comparison only**: CLR(const) vs CLR(std) - no TSS, no coefficient.sm.
3. **Old robustness formula**: `both_sig / (both_sig + const_only)` - two categories only.
4. **SD ribbons**: Used ± 1 SD instead of ± 1 SE.

## Folder structure

- `scripts/` - R script
- `results/` - latest runs per gamma (the ones to look at)
- `old_runs/` - earlier/superseded runs of the same gammas

## Replacement

Use `coeff_sm_reproducibility_spaghetti_cent_mage/` instead.

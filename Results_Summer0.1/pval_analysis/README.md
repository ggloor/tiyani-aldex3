# [DEPRECATED] Old CLR/TSS Constant Regression v1 (cent vs mage, eld vs kin)

**Status: SUPERSEDED** - this is the earliest version (v1) of the constant regression analysis. Superseded by v2.0 (`old_clr_tss_const_regression_v2_all_pairs/`), which was itself superseded by `coeff_sm_regression_analysis_all_pairs/` using the correct `coefficient.sm` from ALDEx3 v1.2.0.

## What was wrong

Same bug as all old constant analyses: custom `clr_const` and `tss_const` functions applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope). Only the slope should be perturbed.

## Differences from later versions

- Single gamma (0.001) - no multi-gamma sweep
- No gamma subfolder nesting (data/plots directly under CLR_constant/ and TSS_constant/)
- v1 script before scatter and effect cross plots were added

## Script

`summary_focused_analysis.R` (copy included; source in `code/`)

## Method (old, incorrect)

- Custom `clr_const` and `tss_const` functions
- Gamma: 0.001 only
- Constants: -1 to 1 by 0.1
- Comparisons: cent vs mage, eld vs kin
- Dataset: Bian et al. 2017 (1117 OTUs)

## Output structure

- `CLR_constant/{data,plots}/` - CLR(const) vs TSS(std)
- `TSS_constant/{data,plots}/` - CLR(std) vs TSS(const)

## Replacement

Use `coeff_sm_regression_analysis_all_pairs/` instead.

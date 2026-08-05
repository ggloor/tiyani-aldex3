# Results_Summer0.1

ALDEx3 differential abundance analysis results on the Bian et al. 2017 gut microbiome dataset (1117 OTUs, 7 age cohorts: cent, eld, kin, mage, mid, pup, you).

## Current analyses (use these)

These use `coefficient.sm` from ALDEx3 v1.2.0 with the correct API: `c.mu = c(0, k)` and `c.cor = diag(c(0, gamma^2))`, which perturbs only the slope coefficient while keeping the intercept at 0.

| Folder | Description | Comparison | Script |
|--------|-------------|------------|--------|
| `coeff_sm_effect_sweep_cent_mage/` | Effect space sweep across constants (-1 to 1) and gammas (0.0001, 0.3, 0.5). 3-panel effect cross + 2-panel scatter per constant. | cent vs mage | `coefficient_sm_effect_sweep.R` |
| `coeff_sm_regression_analysis_all_pairs/` | Regression analysis on red dots (features sig in all three methods). Slope/intercept/R² across constants, gammas 0.1/0.3/0.5. | cent vs mage, eld vs kin | `summary_focused_analysis_3.0.R` |
| `coeff_sm_reproducibility_spaghetti_cent_mage/` | Monte Carlo reproducibility: 24 replicates, spaghetti plots with SE ribbons for robustness, median/mean -log10(padj), agreement counts. | cent vs mage | `clr_vs_clr_reproducibility_3.0.R` |
| `coeff_sm_single_run_test_cent_mage/` | Single-run test (no replicates). Quick verification of the 3-panel effect cross and 5-category system. k=0.3, gamma=0.3. | cent vs mage | `coefficient_sm_TEST.R` |

### Shared conventions across current analyses

- **5 unified categories**: all3 (red), coeff only (orange), CLR only (dodgerblue), TSS only (forestgreen), non-sig (grey75)
- **Robustness formula**: `all_three / (all_three + coeff + clr + tss)` - 1 = all agree, 0 = max disagreement
- **Zero p-value handling**: `min(nonzero)/10` before -log10 transform (Dr. Gloor convention)
- **SE ribbons** (spaghetti plots): `SD / sqrt(n_replicates)`

## ALDEx3 vs edgeR comparison

| Folder | Description |
|--------|-------------|
| `aldex_vs_edger/` | Systematic comparison of ALDEx3 (CLR, TSS) vs edgeR (TMM, TMMwsp, RLE) across all 21 pairwise comparisons. Contains effect cross plots, effect space plots, p-value scatter plots, and Dr. Gloor's preferred single-panel style. See subfolder READMEs. |

Master data-generation scripts for this analysis are at the root level:
- `Updated_Results0.1.R` - runs ALDEx3 + edgeR and saves Rda files
- `Plot_Pval_Scatter_AldexVsEdgeR.R` - p-value scatter plotting (earlier version, gamma=0.3)

## Deprecated analyses (old custom constant functions)

These used custom `clr_const` / `tss_const` functions that incorrectly applied mean k and SD gamma to ALL design matrix coefficients (intercept + slope) via `rnorm(P * ns, k, gamma)`. The intercept should remain at 0 with no uncertainty. All have been superseded by the `coeff_sm_*` folders above.

| Folder | Description | Version |
|--------|-------------|---------|
| `old_clr_tss_const_regression_v2_all_pairs/` | Latest (v2.0) old regression analysis, CLR + TSS constant sides, gammas 0.1/0.3/0.5, cent vs mage + eld vs kin | v2.0 (latest old) |
| `old_clr_const_regression_analysis_all_pairs/` | CLR constant side only, gammas 0.1/0.3/0.5, both pairs | v2.0 (split) |
| `old_tss_const_regression_analysis_all_pairs/` | TSS constant side only, gammas 0.1/0.3/0.5, both pairs | v2.0 (split) |
| `old_clr_tss_const_test_run_cent_mage/` | Test run, gamma 0.3 only, cent vs mage only | v2.0 test |
| `pval_analysis/` | Earliest regression version, gamma 0.001, no multi-gamma sweep | v1 |
| `old_clr_const_gamma_sweep_all_pairs/` | Scatter + effect cross sweep across extreme gammas (1e-5 to 3), all 21 pairs | exploratory |
| `old_clr_const_reproducibility_spaghetti_cent_mage/` | Old spaghetti plots, CLR(const) vs CLR(std) only, SD ribbons | v1 |

## Source scripts

All source scripts are in `code/`. Each result folder contains a copy of its script for reproducibility. See individual folder READMEs for script-to-results mapping.

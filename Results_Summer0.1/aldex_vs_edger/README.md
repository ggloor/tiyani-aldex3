# ALDEx3 vs edgeR Comparison

Systematic comparison of ALDEx3 (CLR and TSS scale models) against edgeR (TMM, TMMwsp, RLE normalizations) across all 21 pairwise age-cohort comparisons from the Bian et al. 2017 gut dataset (1117 OTUs, 7 cohorts).

## Master data script

`Updated_Results0.1.R` (located in `Results_Summer0.1/`) generates all Rda data files stored in `effect_cross_edgeR_allnorms/`. Every plotting script in this folder depends on those Rda files.

- ALDEx3: nsample = 128, gamma = 0.3, CLR and TSS scale models
- edgeR: TMM, TMMwsp, RLE normalizations; both filtered (filterByExpr) and unfiltered runs stored
- FDR threshold: 0.05

## Subfolders

- `effect_cross_edgeR_allnorms/` - Rda data files + unfiltered effect cross plots
- `effect_cross_edgeR_filtered_coords/` - effect cross plots using filtered edgeR coordinates
- `effect_space_aldex_vs_edgeR/` - effect cross plots with both panels in estimate-vs-SE space
- `effect_plot_prof_style/` - Dr. Gloor's preferred single-panel layered effect plot style
- `pval_scatter_aldex_vs_edgeR/` - -log10(FDR) scatter plots (gamma = 0.3)
- `pval_scatter_aldex_vs_edgeR(gamma0)/` - -log10(FDR) scatter plots (gamma = 0)

See individual subfolder READMEs for details.

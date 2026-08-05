# P-value Scatter Plots - ALDEx3 vs edgeR (gamma = 0.3)

-log10(FDR) scatter plots comparing ALDEx3 (CLR or TSS) vs edgeR (TMM, TMMwsp, RLE). Uses filtered OTUs only (those passing filterByExpr). Zero p-values handled with min(nonzero)/10.

- x-axis: -log10(FDR) edgeR
- y-axis: -log10(FDR) ALDEx3
- Red = both sig, Blue = edgeR only, Orange = ALDEx only, Grey = neither

## Script

`Plot_Pval_Scatter_AldexVsEdgeR.R` (copy in this folder)

Depends on Rda files from `effect_cross_edgeR_allnorms/`.

## Plots & Data

- `plots/`: 6 PDFs - CLR/TSS × TMM/TMMwsp/RLE, 21 pages each
- `pval_scatter_summary.csv`: summary statistics per comparison per norm

## Config

gamma = 0.3, FDR < 0.05, all 21 pairwise comparisons

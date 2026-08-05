# P-value Scatter Plots - ALDEx3 vs edgeR (gamma = 0)

Same analysis as `pval_scatter_aldex_vs_edgeR/` but run with gamma = 0 (no scale uncertainty). This tests how ALDEx3 and edgeR compare when ALDEx3 has no gamma perturbation.

No script copy in this folder - produced by the same `Plot_Pval_Scatter_AldexVsEdgeR.R` (see sibling folder) with `gamma_val <- 0`.

## Plots & Data

- `plots/`: 6 PDFs - CLR/TSS × TMM/TMMwsp/RLE, 21 pages each
- `pval_scatter_summary.csv`: summary statistics

## Config

gamma = 0, FDR < 0.05, all 21 pairwise comparisons

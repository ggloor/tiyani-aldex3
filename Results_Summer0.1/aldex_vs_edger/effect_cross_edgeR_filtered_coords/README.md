# Effect Cross Plots - ALDEx3 vs edgeR (filtered edgeR coordinates)

Same 2-panel effect cross layout as `effect_cross_edgeR_allnorms/`, but the edgeR MA panel uses filtered coordinates (kept OTUs only from filterByExpr) instead of unfiltered.

- Left panel: ALDEx3 effect space (all 1117 OTUs), edgeR sig overlay
- Right panel: edgeR MA space (kept OTUs only), ALDEx3 sig overlay (subset to kept)

## Script

`Plot_Filtered_Coords.R` (copy in this folder)

Depends on Rda files from `effect_cross_edgeR_allnorms/`.

## Plots

- `plots0.1/`: 6 PDFs - CLR/TSS × TMM/TMMwsp/RLE, 21 pages each

## Config

gamma = 0.3, FDR < 0.05, all 21 pairwise comparisons

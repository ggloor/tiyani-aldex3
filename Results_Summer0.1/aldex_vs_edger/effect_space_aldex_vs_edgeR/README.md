# Effect Space Plots - ALDEx3 vs edgeR (both panels in effect space)

Puts both ALDEx3 and edgeR panels in the same coordinate system (estimate/logFC vs SE), unlike the previous effect cross plots which used MA space (logFC vs logCPM) for the edgeR panel.

- Left panel: ALDEx3 effect space (estimate vs std_error, all 1117 OTUs)
- Right panel: edgeR effect space (logFC vs SE, filtered/kept OTUs only)

edgeR SE is computed as `|logFC| / sqrt(F_stat)`.

## Script

`Plot_Effect_Space_AldexVsEdgeR.R` (copy in this folder)

Depends on Rda files from `effect_cross_edgeR_allnorms/`.

## Plots

- `plots/`: 6 PDFs - CLR/TSS × TMM/TMMwsp/RLE, 21 pages each

## Config

gamma = 0.3, FDR < 0.05, all 21 pairwise comparisons

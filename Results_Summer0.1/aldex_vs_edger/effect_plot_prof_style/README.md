# Effect Plots - Dr. Gloor's Preferred Style

Single-panel ALDEx3 effect plots with layered significance overlay (Dr. Gloor's style):

- x-axis: `std_error * sqrt(238)`
- y-axis: estimate
- Base layer: all OTUs (black open circles)
- Layer 2: edgeR-significant OTUs (orange, filled)
- Layer 3: ALDEx-significant OTUs (red, filled, smaller) - drawn on top so dual-sig appears red

Produces both filtered and unfiltered edgeR significance versions: 2 ALDEx norms × 3 edgeR norms × 2 filter versions = 12 PDFs, 21 pages each.

## Script

`Plot_Effect_Prof_Style.R` (copy in this folder)

Depends on Rda files from `effect_cross_edgeR_allnorms/`.

## Plots

Script only - no plots generated yet (plots/ folder not created).

## Config

gamma = 0.3, FDR < 0.05, sqrt(n) = sqrt(238), all 21 pairwise comparisons

# CLR Pairwise Effect Plots - All Pairs

Runs ALDEx3 CLR normalization on all 21 pairwise comparisons with two gamma
values and overlays significance on effect cross plots.

## Script

Looped ALDEx 3 Code.R

## What it does

1. Loads all 7 cohorts (cent, eld, kin, mage, mid, pup, you)
2. For each of the 21 pairs, runs ALDEx3 CLR twice:
   - Without gamma: gamma ~ 0 (1e-3)
   - With gamma: gamma = 0.3
3. Plots effect cross (std error vs estimate) per pair with 3 significance layers:
   - Dark green: all features (background)
   - Orange: significant in the no-gamma result (padj < 0.05)
   - Red: significant in the gamma=0.3 result (padj < 0.05)
4. All 21 plots saved to a single multi-page PDF

## Parameters

- Normalization: CLR (clr.sm)
- Gamma values: 1e-3 (no gamma) and 0.3
- nsample: 128
- Significance: padj < 0.05

## Results

- pairwise_analysis_plots_with_clr.pdf - 21-page PDF, one effect cross plot per pair

## Notes

- Orange points show features significant under near-zero gamma (standard CLR)
- Red points show features significant under gamma=0.3
- Comparison of orange vs red shows how gamma affects significance calls

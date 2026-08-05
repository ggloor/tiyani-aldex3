# Effect Cross Plots - ALDEx3 vs edgeR (all norms, unfiltered coords)

Contains the master Rda data files for all 21 pairwise comparisons, plus unfiltered effect cross plots.

## Script

`Updated_Results0.1.R` (in `Results_Summer0.1/`)

## Data

- 21 `.Rda` files (one per pairwise comparison), each containing:
  - `clr`: ALDEx3 CLR results (estimate, std_error, pval_adj)
  - `tss`: ALDEx3 TSS results
  - `edger$TMM`, `edger$TMMwsp`, `edger$RLE`: edgeR results with coords_unfilt, coords_filt, fdr, fdr_unfilt, keep mask

## Plots

- `plots_unfiltered/`: 2-panel effect cross PDFs (ALDEx3 effect space vs edgeR MA space) using unfiltered edgeR coordinates
- 6 PDFs: CLR/TSS × TMM/TMMwsp/RLE, 21 pages each

## Config

nsample = 128, gamma = 0.3, FDR < 0.05, all 7 cohorts (21 pairs)

# Pairwise P-values Data

CSV files containing per-feature p-values, estimates, and std errors from
all 21 pairwise comparisons under different normalizations and gamma values.

## Results

- pairwise_pvalues.csv (Mar 27) - default pairwise p-values
- pairwise_pvalues_clr.csv (Mar 24) - CLR normalization
- pairwise_pvalues_tss.csv (Mar 24) - TSS normalization
- pairwise_pvalues_clr_gamma1e-5.csv - CLR with gamma ~ 0
- pairwise_pvalues_tss_gamma1e-5.csv - TSS with gamma ~ 0
- pairwise_pvalues_with_tss.csv (Feb 24) - early version with TSS
- tss_pairwise_pvalues.csv (Mar 27) - TSS p-values

## Columns (typical)

taxon, estimate, std.error, p.val, p.val.adj, comparison

## Usage

These CSVs are used as input by downstream scripts (e.g., Histograms_Code.R
reads pairwise_pvalues_clr.csv to plot effect size distributions).

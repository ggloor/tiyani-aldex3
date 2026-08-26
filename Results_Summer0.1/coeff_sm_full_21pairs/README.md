# coefficient.sm Full 21-Pair Analysis

Runs coefficient.sm vs CLR(std) vs TSS(std) for all 21 pairwise comparisons
of the 7 Bian et al. age cohorts (cent, eld, kin, mage, mid, pup, you).

## Script

coefficient_sm_full_21pairs.R (copy included in this folder; source in code/)

## Parameters

- Gammas: 0.0001, 0.3, 0.5
- Constants (c.mu[2]): -1 to 1 by 0.1 (21 values)
- nsample: 32
- ALDEx3 API: c.sd (v1.3.0)
  - c.mu = c(0, k) - intercept=0, slope=k
  - c.sd = c(0, gamma) - no intercept uncertainty, gamma on slope

## Cohorts (7 total, 21 pairs)

cent, eld, kin, mage, mid, pup, you

All pairs generated via combn(): cent vs eld, cent vs kin, ..., pup vs you.

## Folder structure

```
coeff_sm_full_21pairs/
  README.md
  coefficient_sm_full_21pairs.R    (script copy)
  gamma_0.0001/
    data/
      summary_stats_gamma0.0001.csv
      slope_table_gamma0.0001.csv
      red_dot_details_gamma0.0001.csv
    plots/
      regression_{pair}_gamma0.0001.pdf   (5 pages: slope coeff-CLR, slope coeff-TSS, R2, robustness, median |effect|)
      scatter_{pair}_gamma0.0001.pdf      (1 page per constant: 2-panel coeff vs CLR + coeff vs TSS)
      effect_{pair}_gamma0.0001.pdf       (1 page per constant: 3-panel effect cross)
  gamma_0.3/
    data/ ...
    plots/ ...
  gamma_0.5/
    data/ ...
    plots/ ...
```

## CSV descriptions

### summary_stats

One row per constant per comparison. Contains:
- robustness (all_three / any_sig)
- category counts (n_all3, n_coeff_only, n_clr_only, n_tss_only, n_none)
- per-method sig counts
- median -log10(padj) per method
- median |effect size| per method
- correlations and MAE between method pairs

### slope_table

Red-dot (all3) regression statistics per constant per comparison.
Direction split: fits separate regressions for dots where coeff is more
significant vs where CLR/TSS is more significant. Reports slope,
intercept, and R-squared for each direction.

### red_dot_details

Per-feature details for all features significant in all three methods.
Includes adjusted p-values, -log10 values, effect estimates, and
standard errors for each method.

## Conventions

- 5 categories: all3 (red), coeff only (orange), CLR only (dodgerblue),
  TSS only (forestgreen), non-sig (grey75)
- Robustness = all_three / (all_three + coeff_only + clr_only + tss_only)
  - 1 = all normalizations agree, 0 = max disagreement
- Zero p-value handling: min(nonzero)/10 before -log10 (Dr. Gloor convention)
- Significance threshold: padj < 0.05

## How to run

```r
source("coefficient_sm_full_21pairs.R")
```

Requires: ALDEx3 (v1.3.0+), data/*.Rda files (7 cohorts).
Runtime: ~several hours depending on system (21 pairs x 21 constants x
3 gammas = 1323 aldex() calls x 3 methods each = 3969 total runs).

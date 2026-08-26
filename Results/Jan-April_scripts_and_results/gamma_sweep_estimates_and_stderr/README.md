# Gamma Sweep - Estimates and Std Error Comparison

Systematic comparison of CLR vs TSS estimates and standard errors across
8 gamma values (1e-5, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5) for all
21 pairwise comparisons.

## Results

### CSVs (per gamma, per normalization)
- clr_gamma_{val}.csv - CLR p-values/estimates at each gamma
- tss_gamma_{val}.csv - TSS p-values/estimates at each gamma

### PDFs - Estimate scatter (CLR vs TSS)
- clr_vs_tss_estimates_gamma_{val}.pdf - scatter of estimates, one per gamma

### PDFs - Std error scatter (CLR vs TSS)
- clr_vs_tss_std-error_gamma_{val}.pdf - scatter of std errors, one per gamma

## Purpose

Shows how the relationship between CLR and TSS changes across gamma values
in both the estimate (effect size) and standard error dimensions. At low
gamma the two methods can diverge substantially; at higher gamma they
tend to converge.

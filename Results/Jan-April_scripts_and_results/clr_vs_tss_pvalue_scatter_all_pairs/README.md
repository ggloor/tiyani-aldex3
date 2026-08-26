# CLR vs TSS P-value Scatter Plots - All Pairs

Scatter plots comparing -log10(padj) between CLR and TSS normalizations
across all 21 pairwise comparisons at different gamma values.

## Results

- clr_vs_tss_(gmma~1e-5)_plots.pdf - gamma ~ 0 (near-zero)
- clr_vs_tss_(gmma=0.1)_plots.pdf - gamma = 0.1
- clr_vs_tss_(gmma=0.3)_plots.pdf - gamma = 0.3
- clr_vs_tss_(gmma=0.4)_plots.pdf - gamma = 0.4
- clr_vs_tss_(gmma=0.5)_plots.pdf - gamma = 0.5
- clr_vs_tss_plots.pdf - default gamma version
- clr_vs_tss_Scatter_plots.pdf - scatter variant
- clr_vs_tss_pvalues_gamma05_TEST.pdf - test run at gamma=0.5
- 0.1clr_vs_tss_(gmma~1e-5)_plots.pdf - variant with 0.1 constant

## Purpose

Shows agreement/disagreement between CLR and TSS normalizations.
Points on the diagonal indicate agreement; deviation shows where the
two methods call different features as significant.

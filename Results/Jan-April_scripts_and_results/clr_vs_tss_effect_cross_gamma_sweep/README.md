# CLR vs TSS Effect Cross - Gamma Sweep

Side-by-side effect cross plots (std error vs estimate) comparing CLR and TSS
normalizations across a range of gamma values.

## Results

- effect_cross_clr_tss_gamma_1e-05.pdf - gamma ~ 0
- effect_cross_clr_tss_gamma01.pdf - gamma = 0.1
- effect_cross_clr_tss_gamma02.pdf - gamma = 0.2
- effect_cross_clr_tss_gamma03.pdf - gamma = 0.3
- effect_cross_clr_tss_gamma04.pdf - gamma = 0.4
- effect_cross_clr_tss_gamma05.pdf - gamma = 0.5
- effect_cross_0.2_clr_tss_gamma05.pdf - variant with 0.2 constant at gamma=0.5

## Purpose

Visualizes how increasing gamma affects the shape of the effect cross
(estimate vs std error funnel) for both CLR and TSS. Higher gamma
adds more prior weight, shrinking estimates toward zero and reducing
the number of significant features.

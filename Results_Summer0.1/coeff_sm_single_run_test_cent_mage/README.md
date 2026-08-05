# Coefficient.sm Single-Run Test (cent vs mage)

Single-run (no replicates) test of coefficient.sm vs CLR(std) vs TSS(std) for cent vs mage. Used to verify the 3-panel effect cross plots and 5-category significance system before running full sweeps.

## Script

`coefficient_sm_TEST.R` (copy included; source in `code/`)

## Method

- **coefficient.sm**: `c.mu = c(0, k)`, `c.cor = diag(c(0, gamma^2))` per ALDEx3 v1.2.0
- **CLR**: `clr.sm` with matching gamma
- **TSS**: `tss.sm` with matching gamma
- nsample = 32, k = 0.3, gamma = 0.3
- Dataset: Bian et al. 2017 (1117 OTUs)
- Comparison: cent vs mage

## Config

- k = 0.3, gamma = 0.3 (single constant, no sweep)

## Plots

- `coeff_test_k0.3_gamma0.3.pdf`
  - Page A: 3-panel effect cross (Coeff, CLR, TSS spaces) with 5 unified categories
  - Page B: 2-panel -log10(FDR) scatter (Coeff vs CLR, Coeff vs TSS)

## Categories

- all3 (red), coeff only (orange), CLR only (dodgerblue), TSS only (forestgreen), non-sig (grey75)
- Robustness = all_three / (all_three + coeff + clr + tss)

## Zero p-value handling

min(nonzero)/10 before -log10 transform (Dr. Gloor convention)

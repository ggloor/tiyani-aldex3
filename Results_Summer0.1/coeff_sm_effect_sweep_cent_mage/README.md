# Coefficient.sm Effect Sweep

Sweeps the coefficient.sm scale model constant (k = c.mu[2]) from -1 to 1 by 0.1 for cent vs mage. Compares coefficient.sm vs CLR(std) vs TSS(std) at each constant value.

## Script

`coefficient_sm_effect_sweep.R` (copy included in this folder; source in `code/`)

## Method

- **coefficient.sm**: `c.mu = c(0, k)`, `c.cor = diag(c(0, gamma^2))` per ALDEx3 v1.2.0
- **CLR**: `clr.sm` with matching gamma
- **TSS**: `tss.sm` with matching gamma
- nsample = 32, dataset = Bian et al. 2017 (1117 OTUs)

## Plots per constant

- **Page A**: 3-panel effect cross (Coeff, CLR, TSS) with unified 5-category coloring and per-panel legends
- **Page B**: 2-panel -log10(FDR) scatter (Coeff vs CLR, Coeff vs TSS)

## Categories

- Red = All three agree (coeff & CLR & TSS significant)
- Orange = Coeff only
- Blue = CLR only
- Green = TSS only
- Grey = Non-significant

## Robustness

`all_three / (all_three + coeff_only + clr_only + tss_only)`. 1 = all normalizations agree, 0 = max disagreement.

## Results

- `results/gamma1e-04/` - gamma = 0.0001 (near-zero uncertainty)
- `results/gamma0.3/` - gamma = 0.3
- `results/gamma0.5/` - gamma = 0.5

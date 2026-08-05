# Coefficient.sm Reproducibility Spaghetti Plots (cent vs mage)

Monte Carlo reproducibility test: runs coefficient.sm vs CLR(std) vs TSS(std) across constants (-1 to 1 by 0.1) with 1 control + 24 replicates. Spaghetti plots show replicate variance with mean +/- SE ribbons.

## Script

`clr_vs_clr_reproducibility_3.0.R` (copy included; source in `code/`)

## Method

- **coefficient.sm**: `c.mu = c(0, k)`, `c.cor = diag(c(0, gamma^2))` per ALDEx3 v1.2.0
- **CLR**: `clr.sm` with matching gamma
- **TSS**: `tss.sm` with matching gamma
- nsample = 32, n_replicates = 24, dataset = Bian et al. 2017 (1117 OTUs)
- Comparison: cent vs mage

## Spaghetti plot pages

1. **Robustness** - `all_three / (all_three + coeff + clr + tss)`, 1 = all agree
2. **Median -log10(padj)** - coefficient model across constants
3. **Mean -log10(padj)** - coefficient model across constants
4. **Three-way agreement count** - features significant in all three methods

## Plot elements

- Black bold line = control (run 1)
- Red dashed line = replicate mean
- Blue thin lines = individual replicates (n=24)
- Blue ribbon = mean +/- 1 SE

## Data

- `coeff_vs_clr_tss_replicates_gamma*.csv` - per-replicate per-constant metrics

## Results

- `results/gamma0.5/` - gamma = 0.5

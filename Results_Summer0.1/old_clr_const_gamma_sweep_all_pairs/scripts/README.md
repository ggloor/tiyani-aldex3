# Scripts

All scripts use the deprecated custom `clr_const` / `tss_const` functions.

- **`clr_const_vs_tss_sweep.R`** - CLR(const) vs TSS scatter plots. Gammas: 1e-5, 1, 2, 3. Constant: -3. nsample = 128. All 21 pairs. Output: `pval_scatter_plots/clr_constant_sweep/`

- **`clr_vs_tss_sweep_gamma_constant.R`** - CLR(std) vs TSS(const) scatter plots. Gamma: 0.001. Constants: -1 to 1 by 0.1. nsample = 128. All 21 pairs. Output: `pval_scatter_plots/gamma_0.001/`

- **`effect_cross_clr_const_sweep.R`** - effect cross plots (estimate vs SE). Gamma: 0.001. Constants: -1 to 1 by 0.1. All 21 pairs. Output: `effect_cross_plots/`

- **`clr_vs_tss_TEST_cent_vs_mage.R`** - quick test script, cent vs mage only. Gamma: 1e-5. Single constant = 3. No PDF output (plots to RStudio).

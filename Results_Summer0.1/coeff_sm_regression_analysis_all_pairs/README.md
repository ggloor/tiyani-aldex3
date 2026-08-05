# Summary Focused Analysis (v3.0)

Three-way comparison of coefficient.sm vs CLR(std) vs TSS(std) across constants (-1 to 1 by 0.1) and multiple gamma values. Includes regression analysis of "all three" (red dot) features.

## Script

`summary_focused_analysis_3.0.R` (copy included; source in `code/`)

## Method

- **coefficient.sm**: `c.mu = c(0, k)`, `c.cor = diag(c(0, gamma^2))` per ALDEx3 v1.2.0
- **CLR**: `clr.sm` with matching gamma
- **TSS**: `tss.sm` with matching gamma
- nsample = 32, dataset = Bian et al. 2017 (1117 OTUs)
- Comparisons: cent vs mage, eld vs kin

## Plots per comparison per gamma

- **regression PDF** - slope, intercept, R² of red-dot regression (Coeff vs CLR, Coeff vs TSS), robustness across constants, median |effect size|
- **scatter PDF** - 2-panel -log10(FDR) scatter (Coeff vs CLR, Coeff vs TSS) per constant
- **effect PDF** - 3-panel effect cross (Coeff, CLR, TSS) per constant with unified 5-category coloring

## Data per gamma

- `summary_stats_gamma*.csv` - robustness, category counts, correlations, MAE, median -log10(padj), median |effect size| per constant per comparison
- `slope_table_gamma*.csv` - regression slope/intercept/R² with direction split for red dots
- `red_dot_details_gamma*.csv` - per-feature details for all-three-significant features

## Robustness

`all_three / (all_three + coeff_only + clr_only + tss_only)`. 1 = all agree, 0 = max disagreement.

## Results

- `gamma_0.1/` - gamma = 0.1
- `gamma_0.3/` - gamma = 0.3
- `gamma_0.5/` - gamma = 0.5

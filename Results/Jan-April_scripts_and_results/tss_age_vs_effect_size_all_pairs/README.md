# TSS Age vs Effect Size - All Pairs

Runs ALDEx3 TSS (gamma=0.3) on all 21 pairwise comparisons of the 7 Bian et al.
age cohorts and summarizes the mean and median effect size per age group.

## Script

Age-mean&median.R

## What it does

1. Loads all 7 cohorts (cent, eld, kin, mage, mid, pup, you)
2. Runs ALDEx3 with TSS normalization (gamma=0.3, nsample=128) for each of the 21 pairs
3. Computes mean and median of the effect estimate per pair
4. Aggregates across all comparisons to get a per-age-group average
5. Orders by age: Kindergarten -> Primary -> Middle School -> Youth -> Middle Age -> Elderly -> Centenarian
6. Plots two figures: age group vs mean effect size, and age group vs median effect size

## Parameters

- Normalization: TSS (tss.sm)
- Gamma: 0.3
- nsample: 128
- All 21 pairwise comparisons

## Results

- age_vs_effect_size.pdf - earlier version output (Feb 25)
- tss_mean_effect_by_age.pdf - later version output (Apr 1)

## Age group ordering

kin (Kindergarten) -> pup (Primary) -> mid (Middle School) -> you (Youth) -> mage (Middle Age) -> eld (Elderly) -> cent (Centenarian)

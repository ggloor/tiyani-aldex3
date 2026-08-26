# CLR Mean vs Median Effect Size - All Pairs

Investigates whether the mean-median difference in effect size estimates
correlates with age difference between cohorts.

## Script

Mean_vs_Median.R

## What it does

1. Loads all 7 cohorts and assigns midpoint ages:
   - kin=4.5, pup=10, mid=13.5, you=22, mage=40, eld=70, cent=97
2. Runs ALDEx3 CLR (gamma=0.3, nsample=128) for each of the 21 pairs
3. For each pair, computes:
   - Median of all effect size estimates
   - Mean of all effect size estimates
   - Difference (mean - median)
   - Absolute age difference between the two cohorts
4. Plots mean-median difference vs age difference with linear regression line
5. Each point labeled with the pair name

## Parameters

- Normalization: CLR (clr.sm)
- Gamma: 0.3
- nsample: 128

## Results

- Rplot_Mean_vs_Median.pdf - scatter plot of mean-median difference vs age difference with regression line (saved manually from R plot window)

## Purpose

Tests whether more age-distant comparisons produce more skewed effect size
distributions (larger mean-median gap), which would indicate asymmetric
differential abundance patterns in more distant age groups.

# CLR ALDEx3 Gamma Test - Kin vs Pup

Early test script for learning the ALDEx3 API and comparing CLR results
with and without gamma on a single pair (kin vs pup).

## Script

ALDEx3_Only_With_Gamma-code.R

## What it does

1. Tests ALDEx3 with a dummy 10x11 matrix first (sanity check)
2. Runs kin vs pup with CLR, gamma ~ 0 (1e-3), nsample=128
3. Plots effect cross (std error vs estimate) with significant features in red
4. Runs again with gamma = 0.3
5. Overlays three layers on the gamma=0.3 effect cross:
   - Dark green: all features
   - Orange: significant without gamma (from the 1e-3 run)
   - Red: significant with gamma=0.3

## Parameters

- Normalization: CLR (clr.sm)
- Gamma values: 1e-3 and 0.3
- nsample: 128
- Single pair: kin vs pup

## Results

No saved output files - plots were displayed to screen only (no pdf()/dev.off() calls).

## Notes

- Contains debugging print statements (Test_Block_1 through 4)
- Comments document early issues with ALDEx3 namespace (ALDEx3::aldex() needed)
  and summary function naming (summary.aldex() changed to summary())
- This script evolved into the looped versions that run all 21 pairs

# Noisy Line Analysis Results

## Overview
This analysis explores how measurement noise affects linear regression fits when both X and Y variables contain uncertainty. The study compares standard (Y on X) and reverse (X on Y) regression methods under different noise conditions.

## Noise Scenarios Tested
1. **Case 1**: No X noise (0.0), moderate Y noise (0.2)
2. **Case 2**: Equal X and Y noise (both 0.2)
3. **Case 3**: More X noise (0.4) than Y noise (0.2)

## Analysis Methods
For each case, we performed:
1. **Forward Fit**: Y = AX + B
2. **Reverse Fit**: X = A'Y + B'

## Results by Case

### Case 1: No X noise, moderate Y noise
- Forward fit: A = 0.9871 ± 0.0556
- Reverse fit: A' = 0.8700
- 1/A = 1.0130
- R² = 0.8588

### Case 2: Equal X and Y noise
- Forward fit: A = 0.8634 ± 0.0593
- Reverse fit: A' = 0.9305
- 1/A = 1.1582
- R² = 0.8034

### Case 3: More X noise than Y noise
- Forward fit: A = 0.6831 ± 0.0668
- Reverse fit: A' = 0.9787
- 1/A = 1.4640
- R² = 0.6685

## Key Findings

### Effect of Noise on Fit Quality
1. R² decreases as noise increases
2. Confidence intervals widen with increased noise
3. The difference between forward and reverse fits increases with noise

### Reversibility (A' vs 1/A)
1. Case 1: Small difference (0.143)
2. Case 2: Moderate difference (0.228)
3. Case 3: Large difference (0.485)

### Choice of Fitting Method
1. When X has no noise: Standard Y on X regression is most appropriate
2. When noise is equal: Both fits are equally valid
3. When X has more noise: Reverse fit might be more appropriate

## Visualizations
The script generates three plots showing:
1. Data points with both forward and reverse fits
2. Visual comparison of fit lines
3. Impact of noise on fit quality

## Implications
1. The choice of regression method matters when both variables have uncertainty
2. Higher noise levels lead to:
   - Lower R² values
   - Larger confidence intervals
   - Greater differences between forward and reverse fits
3. The assumption of "perfect X" in standard regression can lead to biased results when X contains significant noise

## Conclusions
1. Measurement uncertainty in both variables significantly impacts regression results
2. The appropriate choice of regression method depends on relative measurement uncertainties
3. The relationship between forward and reverse fits breaks down with increasing noise
4. Understanding noise levels in both variables is crucial for choosing the appropriate analysis method 
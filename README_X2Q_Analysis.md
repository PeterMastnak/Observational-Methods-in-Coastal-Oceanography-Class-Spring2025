# X2 vs Q Analysis Results

## Overview
This analysis examines the relationship between river flow rate (Q) and the distance metric X2 in an estuary system. The analysis fits both power law and logarithmic relationships to understand how river flow influences salinity intrusion.

## Data
- Q (Dayflow): River flow rate in m³/s
- X2 (CDECX2): Distance in kilometers

## Analysis Methods
1. **Power Law Fit**: X2 = AQ^n
2. **Logarithmic Fit**: log10(X2) = a*log10(Q) + b

## Results

### Power Law Fit (X2 = AQ^n)
- A = 145.82 ± 2.19 (95% CI)
- n = -0.12 ± 0.003 (95% CI)
- Standard Error in X2: 5.20 km

### Logarithmic Fit (log10(X2) = a*log10(Q) + b)
- a = -0.121
- b = 2.166

### Comparison
- Power law exponent (n): -0.1195
- Log fit slope (a): -0.1208
- Difference: 0.0013

## Visualizations
The script generates two plots:
1. **Direct Relationship**: Shows X2 vs Q with power law fit
2. **Log-log Relationship**: Shows log10(X2) vs log10(Q) with linear fit

## Key Findings
1. The negative exponent (-0.12) indicates X2 decreases with increasing flow rate
2. Both fitting methods give very consistent results (difference < 0.002)
3. The relationship is moderately strong with clear physical meaning
4. Standard error of 5.2 km provides a measure of prediction uncertainty

## Physical Interpretation
1. Higher river flow pushes saltwater intrusion downstream (smaller X2)
2. The power law relationship suggests diminishing returns - doubling flow rate doesn't halve X2
3. The consistency between power law and log fits validates the relationship's robustness

## Implications
1. Can predict salinity intrusion distance from flow rate
2. Useful for water management decisions
3. Provides quantitative basis for estuary management
4. Helps understand estuarine response to flow changes 
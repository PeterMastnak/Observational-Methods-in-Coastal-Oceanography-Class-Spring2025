# Oceanographic Data Analysis Tools

This repository contains Python scripts for analyzing oceanographic data using the GSW-Python package (TEOS-10). It was developed as part of a graduate course in Observational Methods in Coastal Oceanography.

## Features

- Exploration of seawater density relationships with pressure, temperature, and salinity
- Analysis of water mass mixing and cabbeling effects
- Generation of synthetic oceanographic data
- Processing and visualization of ocean temperature and salinity profiles
- Comparison of different oceanographic measures (Absolute vs. Practical Salinity, Conservative vs. In-situ Temperature)

## Installation

1. Clone this repository:
```bash
git clone [repository-url]
```

2. Create and activate a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows, use: venv\Scripts\activate
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

## Scripts

- `explore_density_relationships.py`: Investigates how pressure, temperature, and salinity affect seawater density
- `gsw_examples.py`: Demonstrates basic usage of the GSW package and creates T-S diagrams
- `process_ocean_data.py`: Processes and visualizes oceanographic data
- `generate_sample_data.py`: Creates synthetic ocean data for testing and learning

## Usage

Each script can be run independently. For example:

```bash
python explore_density_relationships.py
```

This will generate various plots and a summary report in the current directory.

## Output Files

The scripts generate several visualization files:
- Density relationship plots (PNG format)
- Temperature and salinity profiles
- T-S diagrams
- Summary reports (TXT format)

## Requirements

- Python 3.7+
- gsw-python
- numpy
- matplotlib
- netCDF4
- cmocean

See `requirements.txt` for specific version requirements.

## Contributing

This is an educational project. Feel free to fork and modify for your own use. If you find any bugs or have suggestions, please open an issue.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Developed for the Observational Methods in Coastal Oceanography course
- Uses the GSW-Python package based on TEOS-10 

# C_D Uncertainty Analysis for Palau Tidal Flow

This analysis examines the uncertainty in drag coefficient (C_D) measurements for the Palau tidal flow, considering different pressure sensor separations and tidal flow velocities.

## Evolution of the Analysis Approach

### Initial Approach (Oversimplification)
Initially, we treated the problem using a continuous gradient approach:
- Used direct pressure differences: $C_D = \frac{2\Delta p}{\rho u^2}$
- Assumed continuous pressure gradient
- Did not properly account for the spatial discretization of measurements

This was an oversimplification because:
1. Real measurements are taken at discrete points
2. Ignored the actual spatial arrangement of sensors
3. Did not properly account for the momentum balance

### Correct Approach (Finite Differences)
The proper approach uses finite differences to match how measurements are actually made:
- Uses the finite difference equation: $g\frac{\xi_{i+2}-\xi_i}{x_{i+2}-x_i} = -C_DU_{i+1}|U_{i+1}|/h_{i+1}$
- Explicitly accounts for measurement locations
- Properly represents the spatial discretization
- Matches the actual experimental setup with two pressure sensors

## Methodology

### 1. Finite Difference Implementation
The analysis uses centered finite differences:
- Measurements at positions i and i+2
- Velocity and depth evaluated at i+1
- Properly accounts for the spatial arrangement of measurements

### 2. C_D Calculation
From the finite difference equation:
$C_D = -g h \frac{\xi_{i+2}-\xi_i}{x_{i+2}-x_i} \frac{1}{U_{i+1}|U_{i+1}|}$

Where:
- $\xi_{i+2}, \xi_i$ are water surface elevations at positions i+2 and i
- $x_{i+2}, x_i$ are sensor positions
- $h_{i+1}$ is the water depth at position i+1
- $U_{i+1}$ is the depth-averaged velocity at position i+1

### 3. Uncertainty Analysis
Considers uncertainties in:
- Surface elevation measurements (ξ): ±1mm
- Position measurements (x): ±1cm
- Depth measurements (h): ±1cm
- Velocity measurements (U): ±5cm/s

### 4. Parameter Ranges
- Sensor separations: 0.5 to 5 meters
- Velocities: 0.1 to 1.5 m/s
- Target C_D: 0.025 (midpoint of 0.02-0.03 range)

## Results and Interpretation

### 1. Optimal Sensor Separation
- Recommended separation: approximately 2 meters
- Balances:
  - Need for measurable surface elevation difference
  - Spatial resolution of measurements
  - Practical deployment considerations

### 2. Velocity Effects
- Minimum reliable velocity: 0.5 m/s
- Higher velocities generally yield lower uncertainty
- Very low velocities (< 0.3 m/s) produce unreliable results

### 3. Uncertainty Characteristics
- Velocity uncertainty dominates at low velocities
- Sensor separation becomes more critical at higher velocities
- Surface elevation measurement uncertainty more significant at larger separations

## Recommendations

1. **Measurement Setup**:
   - Use 2-meter sensor separation
   - Ensure accurate position measurements
   - Maintain consistent sensor alignment

2. **Operating Conditions**:
   - Target measurements at velocities > 0.5 m/s
   - Consider multiple measurements during peak tidal flows
   - Avoid very low velocity conditions

3. **Uncertainty Management**:
   - Focus on velocity measurement accuracy
   - Carefully measure sensor positions
   - Consider averaging multiple measurements

## Running the Analysis

1. Required Python packages:
   ```
   numpy
   matplotlib
   scipy
   ```

2. Run the script:
   ```
   python cd_uncertainty_analysis.py
   ```

3. Output:
   - `cd_uncertainty_contour.png`: 2D visualization of uncertainty
   - `cd_uncertainty_vs_separation.png`: Uncertainty vs separation for different velocities 
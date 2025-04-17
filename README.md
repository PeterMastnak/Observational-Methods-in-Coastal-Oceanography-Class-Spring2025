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

## Methodology

The analysis is based on the following approach:

1. **Pressure Difference Measurement**:
   - Using two RBR Quartz pressure sensors
   - Measuring pressure difference (Δp) between sensors
   - Considering sensor separation distances from 0.5m to 5m
   - Analyzing velocities from 0.1 m/s to 1.5 m/s

2. **C_D Calculation**:
   - C_D = 2Δp/(ρu²)
   - Where:
     - Δp is the pressure difference
     - ρ is water density (1025 kg/m³)
     - u is the depth-averaged velocity

3. **Uncertainty Propagation**:
   - Considering uncertainties in:
     - Pressure measurements (Δp)
     - Depth measurements (h)
     - Velocity measurements (u)
   - Using standard error propagation techniques

## Detailed Calculation Method

### 1. Basic C_D Calculation
The drag coefficient (C_D) is calculated using the pressure difference between two sensors:

$C_D = \frac{2\Delta p}{\rho u^2}$

Where:
- Δp is the pressure difference between sensors (Pa)
- ρ is water density (1025 kg/m³)
- u is the depth-averaged velocity (m/s)

### 2. Uncertainty Propagation
We use the standard error propagation formula for uncorrelated variables:

$\sigma_{C_D} = \sqrt{\left(\frac{\partial C_D}{\partial \Delta p}\sigma_{\Delta p}\right)^2 + \left(\frac{\partial C_D}{\partial u}\sigma_u\right)^2}$

The partial derivatives are:
- $\frac{\partial C_D}{\partial \Delta p} = \frac{2}{\rho u^2}$
- $\frac{\partial C_D}{\partial u} = -\frac{4\Delta p}{\rho u^3}$

### 3. Implementation Details
The calculation is implemented in Python with the following key components:

```python
def calculate_cd_uncertainty(delta_p, rho, h, u, delta_h_uncertainty=0.01, delta_p_uncertainty=0.001):
    # Calculate C_D
    cd = 2 * delta_p / (rho * u**2)
    
    # Partial derivatives for uncertainty propagation
    dcd_dp = 2 / (rho * u**2)
    dcd_du = -4 * delta_p / (rho * u**3)
    
    # Calculate total uncertainty using error propagation formula
    cd_uncertainty = np.sqrt(
        (dcd_dp * delta_p_uncertainty)**2 +
        (dcd_du * delta_h_uncertainty)**2
    )
    
    return cd, cd_uncertainty
```

### 4. Analysis Parameters
- Sensor separation range: 0.5 to 5 meters
- Velocity range: 0.1 to 1.5 m/s
- Water density (ρ): 1025 kg/m³
- Assumed uncertainties:
  - Pressure measurement (σ_Δp): 0.001 Pa
  - Depth measurement (σ_h): 0.01 m

### 5. Pressure Difference Calculation
For each combination of separation distance and velocity, we calculate the pressure difference using:

$\Delta p = \frac{1}{2}\rho u^2 \frac{s}{h}$

Where:
- s is the sensor separation distance
- h is the water depth

This relationship comes from the basic hydraulic gradient in the flow.

### 6. Grid Analysis
The script creates a grid of values by:
- Creating arrays of separation distances and velocities
- Computing uncertainties for each combination
- Storing results in a 2D array for plotting

This systematic approach allows us to:
- Identify optimal sensor separation distances
- Understand how uncertainty varies with flow conditions
- Make informed decisions about measurement setup

## Importance of Velocity Uncertainty

The uncertainty in depth-averaged velocity plays a crucial role in estimating CD, and here's why:

### 1. Mathematical Sensitivity
- CD is inversely proportional to velocity squared (u²) in the equation:
  $C_D = \frac{2\Delta p}{\rho u^2}$
- This quadratic relationship means that any uncertainty in velocity gets amplified:
  - A 10% error in velocity measurement leads to approximately 20% error in CD
  - The negative exponent means errors grow larger at lower velocities

### 2. Error Propagation Analysis
- Looking at the partial derivative with respect to velocity:
  $\frac{\partial C_D}{\partial u} = -\frac{4\Delta p}{\rho u^3}$
- The u³ term in the denominator shows that:
  - Velocity uncertainty becomes increasingly critical at lower velocities
  - Small velocity errors have a much larger impact than equivalent pressure measurement errors

### 3. Practical Implications
- At low velocities (< 0.5 m/s):
  - Velocity uncertainty dominates the total uncertainty in CD
  - Measurements become increasingly unreliable
- At higher velocities (> 1.0 m/s):
  - The impact of velocity uncertainty decreases
  - Measurements become more reliable

This analysis suggests that:
1. Accurate velocity measurements are crucial for reliable CD estimates
2. Measurements should ideally be taken during periods of stronger flow
3. Extra care should be taken in velocity measurements at lower flow speeds

## Results and Interpretation

The analysis produces two main visualizations:

1. **Contour Plot** (`cd_uncertainty_contour.png`):
   - Shows C_D uncertainty across different velocities and sensor separations
   - Helps identify optimal sensor separation distances
   - Reveals the relationship between velocity and uncertainty

2. **Uncertainty vs Separation Plot** (`cd_uncertainty_vs_separation.png`):
   - Shows how C_D uncertainty varies with sensor separation
   - Compares different velocity conditions
   - Helps identify optimal separation distances

### Key Findings:

1. **Optimal Sensor Separation**:
   - For typical C_D values (0.02-0.03), a separation distance of approximately 2 meters appears optimal
   - This separation provides a good balance between:
     - Sufficient pressure difference for accurate measurement
     - Minimal uncertainty in the measurements

2. **Velocity Impact**:
   - Higher velocities generally result in lower uncertainty
   - The relationship between velocity and uncertainty is non-linear
   - Very low velocities (< 0.3 m/s) show significantly higher uncertainty

3. **Depth-Averaged Velocity Importance**:
   - Velocity uncertainty has a significant impact on C_D estimation
   - The impact is more pronounced at lower velocities
   - Accurate velocity measurements are crucial for reliable C_D estimates

## Recommendations

1. **Sensor Separation**:
   - Use a separation distance of 2 meters for optimal results
   - This provides good sensitivity while maintaining reasonable uncertainty levels

2. **Measurement Conditions**:
   - Aim to make measurements during higher velocity conditions (> 0.5 m/s)
   - Ensure accurate depth-averaged velocity measurements
   - Consider multiple measurements at different tidal phases

3. **Uncertainty Management**:
   - Focus on reducing velocity measurement uncertainty
   - Use high-precision pressure sensors
   - Consider averaging multiple measurements to reduce random errors

## Running the Analysis

To run the analysis:

1. Ensure you have the required Python packages:
   ```
   numpy
   matplotlib
   scipy
   ```

2. Run the script:
   ```
   python cd_uncertainty_analysis.py
   ```

3. The script will generate two plots:
   - `cd_uncertainty_contour.png`
   - `cd_uncertainty_vs_separation.png` 
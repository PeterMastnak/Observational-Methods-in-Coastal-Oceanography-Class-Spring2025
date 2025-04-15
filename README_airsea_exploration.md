# Air-Sea Flux Exploration with python-airsea

This project provides a comprehensive exploration of air-sea fluxes using the `python-airsea` package, a translation of the original AIR_SEA toolbox for MATLAB. The exploration covers wind stress, heat fluxes, humidity effects, and density relationships.

## Dependencies

- Python 3.9+
- `numpy`
- `matplotlib`
- `gsw` (TEOS-10 Gibbs SeaWater Oceanographic Package)
- `airsea` (can be installed from GitHub: `pip install git+https://github.com/pyoceans/python-airsea.git`)

## Features

This script provides visualizations and analyses of:

1. **Wind Stress Formulations**: Compares different wind stress formulations (Large & Pond, Smith) across a range of wind speeds.

2. **Heat Fluxes**: Explores how sensible and latent heat fluxes vary with air-sea temperature differences.

3. **Humidity Effects**: Analyzes the effect of relative humidity on latent heat flux and evaporation rates.

4. **Density Comparison**: Compares seawater density (using GSW) and air density across a range of temperatures, and examines the density ratio.

5. **Drag Coefficients**: Shows how drag coefficients vary with wind speed for different formulations.

6. **Wind Speed Profile**: Visualizes the logarithmic wind speed profile with height, demonstrating the log law.

## Generated Visualizations

The script generates several informative plots:

- `wind_stress_comparison.png`: Comparison of wind stress formulations
- `heat_flux_comparison.png`: Heat fluxes versus air-sea temperature difference
- `humidity_effect.png`: Effect of humidity on latent heat flux and evaporation
- `density_comparison.png`: Comparison of seawater and air density with temperature
- `drag_coefficient_comparison.png`: Drag coefficients versus wind speed
- `wind_profile.png`: Wind speed profile with height (logarithmic law)

## Usage

Simply run the script:

```bash
python airsea_exploration.py
```

All plots will be generated and saved in the current directory.

## Physical Background

### Wind Stress

Wind stress is the drag force per unit area exerted by the wind on the ocean surface. It is calculated using:

τ = ρₐ × Cd × U²

Where:
- τ is the wind stress (N/m²)
- ρₐ is air density (kg/m³)
- Cd is the drag coefficient (dimensionless)
- U is the wind speed (m/s)

### Heat Fluxes

The script calculates two types of heat fluxes:

1. **Sensible Heat Flux**: Direct transfer of heat due to temperature differences
   H = ρₐ × Cp × Ch × U × (Ts - Ta)

2. **Latent Heat Flux**: Heat transfer associated with evaporation
   E = ρₐ × Lv × Ce × U × (qs - qa)

Where:
- H is sensible heat flux (W/m²)
- E is latent heat flux (W/m²)
- Cp is specific heat capacity of air (J/kg·K)
- Ch is heat transfer coefficient (dimensionless)
- Ce is moisture transfer coefficient (dimensionless)
- Lv is latent heat of vaporization (J/kg)
- qs is saturated specific humidity at sea surface
- qa is specific humidity of air

### Sign Convention

The script follows the oceanographic convention for flux directions:
- **Positive flux**: from ocean to atmosphere (ocean cooling)
- **Negative flux**: from atmosphere to ocean (ocean warming)

## Customization

Users can modify the script to explore different parameter ranges by adjusting the values in the exploration functions. The script provides a foundation for understanding air-sea interaction processes and can be extended for educational or research purposes.

## Acknowledgments

This implementation uses the `python-airsea` package which is a translation of the original AIR_SEA toolbox for MATLAB, and the GSW package for accurate seawater properties calculations based on TEOS-10. 
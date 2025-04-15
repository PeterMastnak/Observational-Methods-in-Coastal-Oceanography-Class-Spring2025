# Air-Sea Flux Toolbox Implementation

This project provides a Python implementation of the air-sea flux toolbox for calculating surface fluxes of momentum, heat, and moisture at the ocean-atmosphere interface. The implementation is based on the COARE 3.0 algorithm with modifications. This project was developed for the Observational Methods in Coastal Oceanography course.

## Files

- **`air_sea_flux_monterey.py`**: Main implementation of air-sea flux calculations with application to Monterey Bay
- **`air_sea_exploration.py`**: Script to explore how different parameters affect air-sea fluxes
- **Output files**:
  - `monterey_bay_fluxes.png`: Visualization of fluxes for Monterey Bay
  - `wind_speed_effect.png`: Analysis of how wind speed affects fluxes
  - `temperature_effect.png`: Analysis of how air-sea temperature difference affects fluxes
  - `humidity_effect.png`: Analysis of how relative humidity affects latent heat flux

## Usage

To calculate air-sea fluxes for Monterey Bay conditions:
```bash
python air_sea_flux_monterey.py
```

To explore the effects of varying different parameters:
```bash
python air_sea_exploration.py
```

## Physical Background

### Surface Fluxes

The air-sea interaction involves the exchange of momentum, heat, and moisture across the air-sea interface. These fluxes are critical for understanding and modeling ocean circulation, climate, and weather patterns.

1. **Momentum Flux (Wind Stress)**: The transfer of momentum from the atmosphere to the ocean, primarily through wind. This drives ocean currents and waves.

2. **Sensible Heat Flux**: The direct heat exchange due to temperature differences between the air and sea surface.

3. **Latent Heat Flux**: The heat exchange associated with phase changes (evaporation/condensation) at the sea surface.

### Sign Convention

In this implementation, we use the oceanographic convention for flux directions:
- **Positive flux**: transfer from ocean to atmosphere (ocean cooling)
- **Negative flux**: transfer from atmosphere to ocean (ocean warming)

### Bulk Aerodynamic Formulas

The bulk aerodynamic formulas are used to calculate surface fluxes:

- **Wind stress**: τ = ρₐ × Cd × U² (where ρₐ is air density, Cd is drag coefficient, U is wind speed)
- **Sensible heat flux**: Hs = ρₐ × Cp × Ch × U × (SST - Tₐ) (where Cp is specific heat, Ch is heat transfer coefficient)
- **Latent heat flux**: Hl = ρₐ × Le × Ce × U × (qs - q) (where Le is latent heat of vaporization, Ce is moisture transfer coefficient)

The transfer coefficients (Cd, Ch, Ce) are not constants but vary with atmospheric stability, wind speed, and surface roughness. They are calculated iteratively using Monin-Obukhov similarity theory.

### Monin-Obukhov Similarity Theory

This theory describes how variables like temperature, humidity, and wind speed vary with height in the atmospheric boundary layer. Key components include:

1. **Monin-Obukhov Length (L)**: A stability parameter that represents the height at which buoyancy forces dominate over mechanical (wind shear) forces.

2. **Stability Functions**: Functions that adjust the transfer coefficients based on atmospheric stability:
   - Stable conditions (L > 0): Air is warmer than the sea surface
   - Unstable conditions (L < 0): Sea surface is warmer than the air

3. **Roughness Lengths**: Parameters that characterize the surface roughness for momentum (z₀), heat (z₀t), and moisture (z₀q).

## Implementation Details

### Key Functions

- **`qsat(T, P)`**: Calculates saturation specific humidity
- **`air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa)`**: Main function that calculates all fluxes
- **`calculate_z0(usr, Ta)`**: Calculates roughness length using Charnock relation
- **`calc_L(usr, tsr, qsr, Ta, theta_v)`**: Calculates Monin-Obukhov length

### Algorithm

1. **Initialization**: Start with neutral stability approximations
2. **Iteration**: 
   - Calculate stability parameters
   - Update transfer coefficients and roughness lengths
   - Recalculate fluxes
   - Update Monin-Obukhov length
3. **Convergence**: Iterate until stability parameters converge

## Example: Monterey Bay

The script provides an example calculation for typical summer afternoon conditions in Monterey Bay:
- Wind speed: 5 m/s
- Sea surface temperature: 15.5°C
- Air temperature: 17.0°C
- Relative humidity: 80%
- Atmospheric pressure: 1015 mb

These conditions typically result in:
- Negative sensible heat flux (heat transfer from atmosphere to ocean)
- Positive latent heat flux (moisture transfer from ocean to atmosphere)
- Stable atmospheric conditions (since air is warmer than sea surface)

## Exploration Results

### Wind Speed Effect
Wind speed significantly affects all fluxes, with:
- Wind stress increasing roughly with the square of wind speed
- Heat and moisture fluxes increasing nearly linearly with wind speed

### Temperature Effect
Air-sea temperature difference strongly affects sensible heat flux and atmospheric stability:
- When air is warmer than sea: negative sensible heat flux, stable atmosphere
- When sea is warmer than air: positive sensible heat flux, unstable atmosphere

### Humidity Effect
Relative humidity primarily affects latent heat flux:
- Higher humidity reduces evaporation and latent heat flux
- Lower humidity increases evaporation and latent heat flux

## References

1. Fairall, C.W., Bradley, E.F., Hare, J.E., Grachev, A.A., Edson, J.B., 2003. Bulk Parameterization of Air–Sea Fluxes: Updates and Verification for the COARE Algorithm. Journal of Climate, 16(4), 571-591.

2. Smith, S.D., 1988. Coefficients for sea surface wind stress, heat flux, and wind profiles as a function of wind speed and temperature. Journal of Geophysical Research, 93(C12), 15467-15472.

3. Pawlowicz, R., et al. Python Sea Water Package (derived from original MATLAB Sea Water package).

## Acknowledgments

This implementation is inspired by various air-sea toolboxes, including the original MATLAB Air-Sea Toolbox, the COARE algorithm, and the Python-AirSea package. The code has been adapted and simplified for educational purposes in the Observational Methods in Coastal Oceanography course. 
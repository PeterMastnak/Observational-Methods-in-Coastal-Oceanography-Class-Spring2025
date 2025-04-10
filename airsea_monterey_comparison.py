#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Exploration of air-sea fluxes in Monterey Bay using the airsea package.
This script calculates fluxes for typical summer afternoon conditions in Monterey Bay
and compares them with the custom implementation from air_sea_flux_monterey.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from airsea import atmosphere as atm
import gsw  # For seawater density calculations

# Set plotting style
plt.style.use('ggplot')

def calc_wind_stress(u, drag_method='largepond'):
    """
    Calculate wind stress using different drag formulations.
    
    Parameters:
    -----------
    u : float or array
        Wind speed (m/s)
    drag_method : str, optional
        Drag coefficient formulation to use:
        'largepond' (Large & Pond, 1981) or 'smith' (Smith, 1988)
        
    Returns:
    --------
    tau : float or array
        Wind stress (N/m²)
    """
    # Air density at standard conditions
    rho_air = 1.22  # kg/m³
    
    # Calculate drag coefficient based on formulation
    if drag_method == 'largepond':
        # Large and Pond (1981)
        if u < 10:
            cd = 1.15e-3
        else:
            cd = 0.49e-3 + 0.065e-3 * u
    elif drag_method == 'smith':
        # Smith (1988)
        cd = (0.61 + 0.063 * u) / 1000
    else:
        raise ValueError(f"Unknown drag method: {drag_method}")
    
    # Calculate wind stress
    tau = rho_air * cd * u**2
    
    return tau, cd

def calc_sensible_heat(u, Ts, Ta, RH=80.0, Pa=1013.0):
    """
    Calculate sensible heat flux using the bulk formula.
    
    Parameters:
    -----------
    u : float or array
        Wind speed (m/s)
    Ts : float or array
        Sea surface temperature (°C)
    Ta : float or array
        Air temperature (°C)
    RH : float or array, optional
        Relative humidity (%), default is 80.0
    Pa : float or array, optional
        Atmospheric pressure (mb), default is 1013.0
        
    Returns:
    --------
    float or array
        Sensible heat flux (W/m²)
    """
    # Air density
    rho_a = atm.air_dens(Ta, RH, Pa)
    
    # Specific heat capacity of air
    Cp = 1004.0  # J/(kg·K)
    
    # Stanton number (heat transfer coefficient)
    # Using a simple approximation
    Ch = 1.3e-3
    
    # Calculate sensible heat flux - positive upward (from sea to air)
    H = rho_a * Cp * Ch * u * (Ts - Ta)
    
    return H

def calc_latent_heat(u, Ts, Ta, RH, Pa=1013.0):
    """
    Calculate latent heat flux using the bulk formula.
    
    Parameters:
    -----------
    u : float or array
        Wind speed (m/s)
    Ts : float or array
        Sea surface temperature (°C)
    Ta : float or array
        Air temperature (°C)
    RH : float or array
        Relative humidity (%)
    Pa : float or array, optional
        Atmospheric pressure (mb), default is 1013.0
        
    Returns:
    --------
    float or array
        Latent heat flux (W/m²)
    """
    # Air density
    rho_a = atm.air_dens(Ta, RH, Pa)
    
    # Calculate specific humidity at saturation (sea surface)
    q_s = atm.qsat(Ts, Pa)
    
    # Calculate specific humidity in air
    q_a = atm.qsat(Ta, Pa) * (RH / 100.0)
    
    # Dalton number (moisture transfer coefficient)
    # Using a simple approximation
    Ce = 1.5e-3
    
    # Latent heat of vaporization
    Lv = 2.5e6  # J/kg
    
    # Calculate latent heat flux - positive upward (from sea to air)
    E = rho_a * Lv * Ce * u * (q_s - q_a)
    
    return E

def monterey_bay_fluxes():
    """Calculate air-sea fluxes for typical summer afternoon conditions in Monterey Bay."""
    
    # Input values for a summer afternoon in Monterey Bay
    # Based on NOAA buoy data for Monterey Bay in summer
    u = 5.0        # Wind speed (m/s) - moderate afternoon sea breeze
    SST = 15.5     # Sea surface temperature (°C) - typical summer value
    Ta = 17.0      # Air temperature (°C) - warmer than sea surface due to land warming
    RH = 80.0      # Relative humidity (%) - typical marine environment
    Pa = 1015.0    # Atmospheric pressure (mb)
    
    # Calculate specific humidity at saturation (sea surface) and in air
    q_s = atm.qsat(SST, Pa)
    q_a = atm.qsat(Ta, Pa) * (RH / 100.0)
    
    # Calculate wind stress using custom function
    tau, cd = calc_wind_stress(u, drag_method='largepond')
    
    # Calculate heat fluxes
    Hs = calc_sensible_heat(u, SST, Ta, RH, Pa)
    Hl = calc_latent_heat(u, SST, Ta, RH, Pa)
    
    # Calculate evaporation rate (mm/day) from latent heat flux
    Lv = 2.5e6  # J/kg, latent heat of vaporization
    evap_rate = max(0, Hl) / Lv  # kg/m²/s (only consider positive fluxes for evaporation)
    evap_mm_per_day = evap_rate * 86400  # mm/day
    
    # Print input parameters and results
    print("\nMonterey Bay Summer Afternoon Conditions:")
    print("-" * 50)
    print(f"Input Parameters:")
    print(f"  Wind speed: {u:.1f} m/s")
    print(f"  Sea surface temperature: {SST:.1f} °C")
    print(f"  Air temperature: {Ta:.1f} °C")
    print(f"  Relative humidity: {RH:.1f} %")
    print(f"  Specific humidity at sea surface: {q_s*1000:.2f} g/kg")
    print(f"  Specific humidity in air: {q_a*1000:.2f} g/kg")
    print(f"  Atmospheric pressure: {Pa:.1f} mb")
    print()
    
    print("Calculated Fluxes (using python-airsea):")
    print(f"  Wind stress (τ): {tau*1000:.2f} mN/m²")
    print(f"  Drag coefficient (Cd): {cd*1000:.2f} × 10⁻³")
    print(f"  Sensible heat flux (Hs): {Hs:.2f} W/m²")
    print(f"  Latent heat flux (Hl): {Hl:.2f} W/m²")
    print(f"  Total turbulent heat flux (Hs + Hl): {Hs + Hl:.2f} W/m²")
    print(f"  Evaporation rate: {evap_mm_per_day:.2f} mm/day")
    
    # Interpretation of fluxes
    if Hs > 0:
        print(f"  -> Sensible heat is transferring from ocean to atmosphere (ocean cooling)")
    else:
        print(f"  -> Sensible heat is transferring from atmosphere to ocean (ocean warming)")
        
    if Hl > 0:
        print(f"  -> Latent heat is transferring from ocean to atmosphere (evaporation)")
    else:
        print(f"  -> Latent heat is transferring from atmosphere to ocean (condensation)")
    
    # Return the results as a dictionary for further analysis or plotting
    return {
        'u': u,
        'SST': SST,
        'Ta': Ta,
        'RH': RH,
        'Pa': Pa,
        'tau': tau,
        'Cd': cd,
        'Hs': Hs,
        'Hl': Hl,
        'total_flux': Hs + Hl,
        'evap_rate': evap_mm_per_day
    }

def plot_monterey_fluxes(fluxes):
    """Create plots of Monterey Bay air-sea fluxes."""
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # Heat fluxes plot
    fluxes_names = ['Sensible Heat', 'Latent Heat']
    fluxes_values = [fluxes['Hs'], fluxes['Hl']]
    
    # Add total as third bar
    fluxes_names.append('Total Heat')
    fluxes_values.append(fluxes['total_flux'])
    
    # Create bar chart with conditional colors (positive/negative)
    colors = ['blue' if x < 0 else 'red' for x in fluxes_values]
    axs[0].bar(fluxes_names, fluxes_values, color=colors, alpha=0.7)
    axs[0].set_ylabel('Heat Flux (W/m²)')
    axs[0].set_title('Turbulent Heat Fluxes in Monterey Bay\n(positive = ocean to atmosphere)')
    axs[0].axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Add flux values as text
    for i, v in enumerate(fluxes_values):
        axs[0].text(i, v + np.sign(v)*5, f"{v:.1f}", ha='center')
    
    # Wind parameters plot
    param_names = ['Wind Speed (m/s)', 'Wind Stress (mN/m²)', 'Cd (×10⁻³)']
    param_values = [fluxes['u'], fluxes['tau']*1000, fluxes['Cd']*1000]
    
    axs[1].bar(param_names, param_values, color='green', alpha=0.7)
    axs[1].set_title('Wind Parameters in Monterey Bay')
    
    # Add values as text
    for i, v in enumerate(param_values):
        axs[1].text(i, v + 0.05*v, f"{v:.2f}", ha='center')
    
    plt.tight_layout()
    plt.savefig('figures/monterey/monterey_bay_airsea_fluxes.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'figures/monterey/monterey_bay_airsea_fluxes.png'")
    
def explore_monterey_parameters():
    """Explore variations around typical Monterey Bay parameters."""
    
    # Base parameters
    base_u = 5.0        # Wind speed (m/s)
    base_SST = 15.5     # Sea surface temperature (°C)
    base_Ta = 17.0      # Air temperature (°C)
    base_RH = 80.0      # Relative humidity (%)
    base_Pa = 1015.0    # Atmospheric pressure (mb)
    
    print("\nExploring variations in Monterey Bay parameters:")
    print("-" * 50)
    
    # 1. Wind speed variations
    print("\n1. Effect of wind speed variations:")
    wind_speeds = np.linspace(2.0, 10.0, 5)  # 2 to 10 m/s
    hs_values_wind = []
    hl_values_wind = []
    
    for u in wind_speeds:
        Hs = calc_sensible_heat(u, base_SST, base_Ta, base_RH, base_Pa)
        Hl = calc_latent_heat(u, base_SST, base_Ta, base_RH, base_Pa)
        hs_values_wind.append(Hs)
        hl_values_wind.append(Hl)
        print(f"  Wind speed {u:.1f} m/s: Sensible flux = {Hs:.2f} W/m², Latent flux = {Hl:.2f} W/m²")
    
    # 2. SST variations
    print("\n2. Effect of sea surface temperature variations:")
    sst_values = np.linspace(13.5, 17.5, 5)  # 13.5 to 17.5 °C
    hs_values_sst = []
    hl_values_sst = []
    
    for sst in sst_values:
        Hs = calc_sensible_heat(base_u, sst, base_Ta, base_RH, base_Pa)
        Hl = calc_latent_heat(base_u, sst, base_Ta, base_RH, base_Pa)
        hs_values_sst.append(Hs)
        hl_values_sst.append(Hl)
        print(f"  SST {sst:.1f} °C: Sensible flux = {Hs:.2f} W/m², Latent flux = {Hl:.2f} W/m²")
    
    # 3. Air temperature variations
    print("\n3. Effect of air temperature variations:")
    ta_values = np.linspace(15.0, 19.0, 5)  # 15 to 19 °C
    hs_values_ta = []
    hl_values_ta = []
    
    for ta in ta_values:
        Hs = calc_sensible_heat(base_u, base_SST, ta, base_RH, base_Pa)
        Hl = calc_latent_heat(base_u, base_SST, ta, base_RH, base_Pa)
        hs_values_ta.append(Hs)
        hl_values_ta.append(Hl)
        print(f"  Air temp {ta:.1f} °C: Sensible flux = {Hs:.2f} W/m², Latent flux = {Hl:.2f} W/m²")
    
    # 4. Relative humidity variations
    print("\n4. Effect of relative humidity variations:")
    rh_values = np.linspace(60.0, 100.0, 5)  # 60% to 100%
    hs_values_rh = []
    hl_values_rh = []
    
    for rh in rh_values:
        Hs = calc_sensible_heat(base_u, base_SST, base_Ta, rh, base_Pa)
        Hl = calc_latent_heat(base_u, base_SST, base_Ta, rh, base_Pa)
        hs_values_rh.append(Hs)
        hl_values_rh.append(Hl)
        print(f"  RH {rh:.1f} %: Sensible flux = {Hs:.2f} W/m², Latent flux = {Hl:.2f} W/m²")
    
    # Create multi-panel plot of parameter variations
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Wind speed variations
    axs[0, 0].plot(wind_speeds, hs_values_wind, 'r-', linewidth=2, label='Sensible Heat Flux')
    axs[0, 0].plot(wind_speeds, hl_values_wind, 'b-', linewidth=2, label='Latent Heat Flux')
    axs[0, 0].set_xlabel('Wind Speed (m/s)')
    axs[0, 0].set_ylabel('Heat Flux (W/m²)')
    axs[0, 0].set_title('Effect of Wind Speed')
    axs[0, 0].legend()
    axs[0, 0].grid(True, alpha=0.3)
    
    # SST variations
    axs[0, 1].plot(sst_values, hs_values_sst, 'r-', linewidth=2, label='Sensible Heat Flux')
    axs[0, 1].plot(sst_values, hl_values_sst, 'b-', linewidth=2, label='Latent Heat Flux')
    axs[0, 1].set_xlabel('Sea Surface Temperature (°C)')
    axs[0, 1].set_ylabel('Heat Flux (W/m²)')
    axs[0, 1].set_title('Effect of Sea Surface Temperature')
    axs[0, 1].legend()
    axs[0, 1].grid(True, alpha=0.3)
    
    # Air temperature variations
    axs[1, 0].plot(ta_values, hs_values_ta, 'r-', linewidth=2, label='Sensible Heat Flux')
    axs[1, 0].plot(ta_values, hl_values_ta, 'b-', linewidth=2, label='Latent Heat Flux')
    axs[1, 0].set_xlabel('Air Temperature (°C)')
    axs[1, 0].set_ylabel('Heat Flux (W/m²)')
    axs[1, 0].set_title('Effect of Air Temperature')
    axs[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)  # Zero line
    axs[1, 0].legend()
    axs[1, 0].grid(True, alpha=0.3)
    
    # Relative humidity variations
    axs[1, 1].plot(rh_values, hs_values_rh, 'r-', linewidth=2, label='Sensible Heat Flux')
    axs[1, 1].plot(rh_values, hl_values_rh, 'b-', linewidth=2, label='Latent Heat Flux')
    axs[1, 1].set_xlabel('Relative Humidity (%)')
    axs[1, 1].set_ylabel('Heat Flux (W/m²)')
    axs[1, 1].set_title('Effect of Relative Humidity')
    axs[1, 1].legend()
    axs[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/monterey/monterey_bay_parameter_variations.png', dpi=300, bbox_inches='tight')
    print("\nPlot saved as 'figures/monterey/monterey_bay_parameter_variations.png'")
    
def main():
    """Main function to run the Monterey Bay air-sea flux analysis."""
    print("\nAir-Sea Flux Analysis for Monterey Bay using python-airsea")
    print("=" * 60)
    
    # Calculate fluxes for typical Monterey Bay conditions
    fluxes = monterey_bay_fluxes()
    
    # Plot the results
    plot_monterey_fluxes(fluxes)
    
    # Explore variations in parameters
    explore_monterey_parameters()
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main() 
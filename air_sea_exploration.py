#!/usr/bin/env python3
"""
Air-Sea Toolbox - Parameter Exploration
Investigate how air-sea fluxes change with different input parameters
"""

import numpy as np
import matplotlib.pyplot as plt
from air_sea_flux_monterey import air_sea_fluxes, qsat

# Constants for plotting
RHO_AIR = 1.22  # kg/m^3, density of air at sea level
LE = 2.5e6  # J/kg, latent heat of evaporation

def explore_wind_speed_effect():
    """
    Explore how wind speed affects air-sea fluxes
    """
    # Fixed parameters
    zu = 2.0       # Height of wind measurement (m)
    SST = 15.5     # Sea surface temperature (°C)
    Ta = 17.0      # Air temperature (°C)
    za = 2.0       # Height of temperature measurement (m)
    RH = 80.0      # Relative humidity (%)
    zq = 2.0       # Height of humidity measurement (m)
    Pa = 1015.0    # Atmospheric pressure (mb)
    
    # Calculate specific humidity
    q_sat_air = qsat(Ta, Pa)
    q = RH / 100.0 * q_sat_air
    
    # Wind speeds to explore
    wind_speeds = np.linspace(1.0, 15.0, 15)  # 1 to 15 m/s
    
    # Arrays to store results
    tau_values = []
    hs_values = []
    hl_values = []
    cd_values = []
    
    # Calculate fluxes for each wind speed
    for u in wind_speeds:
        fluxes = air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa)
        tau_values.append(fluxes['tau'] * 1000)  # Convert to mN/m²
        hs_values.append(fluxes['Hs'])
        hl_values.append(fluxes['Hl'])
        cd_values.append(fluxes['Cd'] * 1000)  # Convert to 10^-3
    
    # Create multi-panel figure
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot wind stress vs wind speed
    axs[0, 0].plot(wind_speeds, tau_values, 'b-', linewidth=2)
    axs[0, 0].set_xlabel('Wind Speed (m/s)')
    axs[0, 0].set_ylabel('Wind Stress (mN/m²)')
    axs[0, 0].set_title('Wind Stress vs Wind Speed')
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot sensible heat flux vs wind speed
    axs[0, 1].plot(wind_speeds, hs_values, 'r-', linewidth=2)
    axs[0, 1].set_xlabel('Wind Speed (m/s)')
    axs[0, 1].set_ylabel('Sensible Heat Flux (W/m²)')
    axs[0, 1].set_title('Sensible Heat Flux vs Wind Speed')
    axs[0, 1].grid(True, alpha=0.3)
    
    # Plot latent heat flux vs wind speed
    axs[1, 0].plot(wind_speeds, hl_values, 'g-', linewidth=2)
    axs[1, 0].set_xlabel('Wind Speed (m/s)')
    axs[1, 0].set_ylabel('Latent Heat Flux (W/m²)')
    axs[1, 0].set_title('Latent Heat Flux vs Wind Speed')
    axs[1, 0].grid(True, alpha=0.3)
    
    # Plot drag coefficient vs wind speed
    axs[1, 1].plot(wind_speeds, cd_values, 'k-', linewidth=2)
    axs[1, 1].set_xlabel('Wind Speed (m/s)')
    axs[1, 1].set_ylabel('Drag Coefficient (×10⁻³)')
    axs[1, 1].set_title('Drag Coefficient vs Wind Speed')
    axs[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/general_airsea/wind_speed_effect.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'figures/general_airsea/wind_speed_effect.png'")

def explore_temperature_effect():
    """
    Explore how air-sea temperature difference affects fluxes
    """
    # Fixed parameters
    u = 5.0        # Wind speed (m/s)
    zu = 2.0       # Height of wind measurement (m)
    SST = 15.5     # Base sea surface temperature (°C)
    za = 2.0       # Height of temperature measurement (m)
    RH = 80.0      # Relative humidity (%)
    zq = 2.0       # Height of humidity measurement (m)
    Pa = 1015.0    # Atmospheric pressure (mb)
    
    # Temperature differences to explore (-5 to +5 °C relative to SST)
    temp_diffs = np.linspace(-5.0, 5.0, 11)
    air_temps = SST + temp_diffs
    
    # Arrays to store results
    hs_values = []
    hl_values = []
    total_heat_values = []
    stability_values = []
    
    # Calculate fluxes for each air temperature
    for Ta in air_temps:
        # Recalculate specific humidity for each air temperature
        q_sat_air = qsat(Ta, Pa)
        q = RH / 100.0 * q_sat_air
        
        fluxes = air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa)
        hs_values.append(fluxes['Hs'])
        hl_values.append(fluxes['Hl'])
        total_heat_values.append(fluxes['Hs'] + fluxes['Hl'])
        stability_values.append(fluxes['L'])
    
    # Create multi-panel figure
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot sensible heat flux vs temperature difference
    axs[0, 0].plot(temp_diffs, hs_values, 'r-', linewidth=2)
    axs[0, 0].set_xlabel('Air-Sea Temperature Difference (°C)')
    axs[0, 0].set_ylabel('Sensible Heat Flux (W/m²)')
    axs[0, 0].set_title('Sensible Heat Flux vs Temperature Difference')
    axs[0, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axs[0, 0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot latent heat flux vs temperature difference
    axs[0, 1].plot(temp_diffs, hl_values, 'g-', linewidth=2)
    axs[0, 1].set_xlabel('Air-Sea Temperature Difference (°C)')
    axs[0, 1].set_ylabel('Latent Heat Flux (W/m²)')
    axs[0, 1].set_title('Latent Heat Flux vs Temperature Difference')
    axs[0, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axs[0, 1].axvline(x=0, color='k', linestyle='--', alpha=0.3)
    axs[0, 1].grid(True, alpha=0.3)
    
    # Plot total heat flux vs temperature difference
    axs[1, 0].plot(temp_diffs, total_heat_values, 'm-', linewidth=2)
    axs[1, 0].set_xlabel('Air-Sea Temperature Difference (°C)')
    axs[1, 0].set_ylabel('Total Heat Flux (W/m²)')
    axs[1, 0].set_title('Total Heat Flux vs Temperature Difference')
    axs[1, 0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axs[1, 0].axvline(x=0, color='k', linestyle='--', alpha=0.3)
    axs[1, 0].grid(True, alpha=0.3)
    
    # Plot Monin-Obukhov length vs temperature difference
    axs[1, 1].plot(temp_diffs, stability_values, 'b-', linewidth=2)
    axs[1, 1].set_xlabel('Air-Sea Temperature Difference (°C)')
    axs[1, 1].set_ylabel('Monin-Obukhov Length (m)')
    axs[1, 1].set_title('Stability vs Temperature Difference')
    axs[1, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axs[1, 1].axvline(x=0, color='k', linestyle='--', alpha=0.3)
    axs[1, 1].grid(True, alpha=0.3)
    
    # Add text annotations for stability regimes
    axs[1, 1].text(2, np.max(stability_values)*0.8, 'Stable', fontsize=12, ha='center')
    axs[1, 1].text(-2, np.min(stability_values)*0.8, 'Unstable', fontsize=12, ha='center')
    
    plt.tight_layout()
    plt.savefig('figures/general_airsea/temperature_effect.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'figures/general_airsea/temperature_effect.png'")

def explore_humidity_effect():
    """
    Explore how relative humidity affects latent heat flux
    """
    # Fixed parameters
    u = 5.0        # Wind speed (m/s)
    zu = 2.0       # Height of wind measurement (m)
    SST = 15.5     # Sea surface temperature (°C)
    Ta = 17.0      # Air temperature (°C)
    za = 2.0       # Height of temperature measurement (m)
    zq = 2.0       # Height of humidity measurement (m)
    Pa = 1015.0    # Atmospheric pressure (mb)
    
    # Relative humidity values to explore
    rh_values = np.linspace(40.0, 100.0, 13)  # 40% to 100%
    
    # Arrays to store results
    hl_values = []
    evap_values = []
    
    # Calculate fluxes for each relative humidity
    for RH in rh_values:
        # Calculate specific humidity
        q_sat_air = qsat(Ta, Pa)
        q = RH / 100.0 * q_sat_air
        
        fluxes = air_sea_fluxes(u, zu, SST, Ta, za, q, zq, Pa)
        hl_values.append(fluxes['Hl'])
        
        # Calculate evaporation rate (mm/day)
        latent_heat_flux = abs(fluxes['Hl']) if fluxes['Hl'] > 0 else 0
        evap_rate = latent_heat_flux / LE  # kg/m²/s
        evap_mm_per_day = evap_rate * 86400  # mm/day
        evap_values.append(evap_mm_per_day)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot latent heat flux vs relative humidity
    ax1.plot(rh_values, hl_values, 'g-', linewidth=2)
    ax1.set_xlabel('Relative Humidity (%)')
    ax1.set_ylabel('Latent Heat Flux (W/m²)')
    ax1.set_title('Latent Heat Flux vs Relative Humidity')
    ax1.grid(True, alpha=0.3)
    
    # Plot evaporation rate vs relative humidity
    ax2.plot(rh_values, evap_values, 'b-', linewidth=2)
    ax2.set_xlabel('Relative Humidity (%)')
    ax2.set_ylabel('Evaporation Rate (mm/day)')
    ax2.set_title('Evaporation Rate vs Relative Humidity')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/general_airsea/humidity_effect.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'figures/general_airsea/humidity_effect.png'")

def main():
    """Main function to run all explorations"""
    print("Exploring Air-Sea Flux Parameter Relationships")
    print("=" * 50)
    
    print("\n1. Exploring effect of wind speed on air-sea fluxes...")
    explore_wind_speed_effect()
    
    print("\n2. Exploring effect of air-sea temperature difference on fluxes...")
    explore_temperature_effect()
    
    print("\n3. Exploring effect of relative humidity on latent heat flux...")
    explore_humidity_effect()
    
    print("\nAll explorations complete!")

if __name__ == "__main__":
    main()

 
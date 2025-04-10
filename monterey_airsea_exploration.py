#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Exploration of air-sea fluxes in Monterey Bay using the airsea package.
This script applies the typical Monterey Bay parameters to the airsea_exploration.py
functionality, focusing on the specific conditions of Monterey Bay.
"""

import numpy as np
import matplotlib.pyplot as plt
from airsea import atmosphere as atm
import gsw  # For comparison with GSW density calculations

# Set plotting style
plt.style.use('ggplot')

# Monterey Bay parameters (typical summer afternoon)
MONTEREY_WIND_SPEED = 5.0        # Wind speed (m/s)
MONTEREY_SST = 15.5              # Sea surface temperature (°C)
MONTEREY_AIR_TEMP = 17.0         # Air temperature (°C)
MONTEREY_RH = 80.0               # Relative humidity (%)
MONTEREY_PRESSURE = 1015.0       # Atmospheric pressure (mb)

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
        if np.isscalar(u):
            # Handle scalar values
            if u < 10:
                cd = 1.15e-3
            else:
                cd = 0.49e-3 + 0.065e-3 * u
        else:
            # Handle arrays
            cd = np.ones_like(u) * 1.15e-3
            mask = u >= 10
            cd[mask] = 0.49e-3 + 0.065e-3 * u[mask]
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

def explore_wind_stress_formulations(base_speed=MONTEREY_WIND_SPEED):
    """Compare different wind stress formulations around Monterey Bay conditions."""
    
    # Create a range of wind speeds to explore (start from 0.1 to avoid divide by zero)
    wind_speeds = np.linspace(0.1, 25, 100)  # 0.1 to 25 m/s
    
    # Calculate wind stress using different formulations
    tau_lp, cd_lp = calc_wind_stress(wind_speeds, drag_method='largepond')
    tau_s, cd_s = calc_wind_stress(wind_speeds, drag_method='smith')
    
    # Mark Monterey Bay conditions
    tau_monterey, cd_monterey = calc_wind_stress(base_speed, drag_method='largepond')
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot wind stress vs wind speed
    ax1.plot(wind_speeds, tau_lp, 'b-', linewidth=2, label='Large & Pond (1981)')
    ax1.plot(wind_speeds, tau_s, 'r-', linewidth=2, label='Smith (1988)')
    ax1.scatter(base_speed, tau_monterey, color='green', s=100, marker='*', 
               label='Monterey Bay condition')
    
    # Add labels and title for first subplot
    ax1.set_xlabel('Wind Speed (m/s)', fontsize=12)
    ax1.set_ylabel('Wind Stress (N/m²)', fontsize=12)
    ax1.set_title('Wind Stress vs. Wind Speed', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Plot drag coefficient vs wind speed
    ax2.plot(wind_speeds, cd_lp * 1000, 'b-', linewidth=2, label='Large & Pond (1981)')
    ax2.plot(wind_speeds, cd_s * 1000, 'r-', linewidth=2, label='Smith (1988)')
    ax2.scatter(base_speed, cd_monterey * 1000, color='green', s=100, marker='*',
               label='Monterey Bay condition')
    
    # Add labels and title for second subplot
    ax2.set_xlabel('Wind Speed (m/s)', fontsize=12)
    ax2.set_ylabel('Drag Coefficient (Cd × 1000)', fontsize=12)
    ax2.set_title('Drag Coefficient vs. Wind Speed', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('figures/general_airsea/figures/monterey/monterey_wind_stress_comparison.png', dpi=300)
    plt.close()
    
    print("Wind stress comparison completed and saved to 'figures/general_airsea/figures/monterey/monterey_wind_stress_comparison.png'")

def explore_heat_fluxes():
    """Explore sensible and latent heat fluxes around Monterey Bay conditions."""
    
    # Create a range of temperature differences
    delta_t = np.linspace(-5, 5, 100)  # From -5°C to +5°C difference
    
    # Fixed parameters based on Monterey Bay
    wind_speed = MONTEREY_WIND_SPEED
    sea_temp = MONTEREY_SST
    base_air_temp = MONTEREY_AIR_TEMP
    pressure = MONTEREY_PRESSURE
    humidity = MONTEREY_RH
    
    # Calculate the air temperatures
    air_temp = sea_temp + delta_t
    
    # Calculate heat fluxes
    sensible_heat = np.array([calc_sensible_heat(wind_speed, sea_temp, t, humidity, pressure) for t in air_temp])
    latent_heat = np.array([calc_latent_heat(wind_speed, sea_temp, t, humidity, pressure) for t in air_temp])
    
    # Total heat flux
    total_heat = sensible_heat + latent_heat
    
    # Mark Monterey Bay condition (delta_t = Ta - SST = 17.0 - 15.5 = 1.5)
    monterey_delta_t = base_air_temp - sea_temp
    monterey_sensible = calc_sensible_heat(wind_speed, sea_temp, base_air_temp, humidity, pressure)
    monterey_latent = calc_latent_heat(wind_speed, sea_temp, base_air_temp, humidity, pressure)
    monterey_total = monterey_sensible + monterey_latent
    
    # Create the plot
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # Plot sensible heat flux
    ax1.plot(delta_t, sensible_heat, 'b-', linewidth=2, label='Sensible Heat Flux')
    # Plot latent heat flux
    ax1.plot(delta_t, latent_heat, 'r-', linewidth=2, label='Latent Heat Flux')
    # Plot total heat flux
    ax1.plot(delta_t, total_heat, 'g-', linewidth=2, label='Total Heat Flux')
    
    # Mark Monterey Bay conditions
    ax1.scatter(monterey_delta_t, monterey_sensible, color='blue', s=100, marker='*')
    ax1.scatter(monterey_delta_t, monterey_latent, color='red', s=100, marker='*')
    ax1.scatter(monterey_delta_t, monterey_total, color='green', s=100, marker='*')
    ax1.axvline(x=monterey_delta_t, color='k', linestyle='--', alpha=0.5, label='Monterey Bay condition')
    
    # Add labels and title
    ax1.set_xlabel('Temperature Difference (Sea - Air) °C', fontsize=12)
    ax1.set_ylabel('Heat Flux (W/m²)', fontsize=12)
    ax1.set_title('Heat Fluxes vs. Air-Sea Temperature Difference (Monterey Bay)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)
    
    # Add a text box explaining sign convention
    textstr = "Positive flux: from ocean to atmosphere (ocean cooling)\nNegative flux: from atmosphere to ocean (ocean warming)"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Add Monterey Bay parameter text box
    monterey_str = f"Monterey Bay Parameters:\nWind speed: {wind_speed} m/s\nSST: {sea_temp}°C\nAir temp: {base_air_temp}°C\nRH: {humidity}%"
    props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.5)
    ax1.text(0.75, 0.25, monterey_str, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('figures/general_airsea/figures/monterey/monterey_heat_flux_comparison.png', dpi=300)
    plt.close()
    
    print("Heat flux comparison completed and saved to 'figures/general_airsea/figures/monterey/monterey_heat_flux_comparison.png'")

def explore_humidity_effect():
    """Explore the effect of relative humidity on latent heat flux around Monterey Bay conditions."""
    
    # Create a range of relative humidity values
    humidity = np.linspace(40, 100, 100)  # From 40% to 100%
    
    # Fixed parameters based on Monterey Bay
    wind_speed = MONTEREY_WIND_SPEED
    sea_temp = MONTEREY_SST
    air_temp = MONTEREY_AIR_TEMP
    pressure = MONTEREY_PRESSURE
    base_rh = MONTEREY_RH
    
    # Calculate latent heat flux
    latent_heat = np.array([calc_latent_heat(wind_speed, sea_temp, air_temp, h, pressure) for h in humidity])
    
    # Calculate evaporation rate (mm/day)
    evaporation = latent_heat / (2.5e6) * 86400  # Convert W/m² to mm/day using latent heat of vaporization
    
    # Mark Monterey Bay condition
    monterey_latent = calc_latent_heat(wind_speed, sea_temp, air_temp, base_rh, pressure)
    monterey_evap = monterey_latent / (2.5e6) * 86400
    
    # Create the plot with two y-axes
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Plot latent heat flux
    ax1.plot(humidity, latent_heat, 'b-', linewidth=2)
    ax1.set_xlabel('Relative Humidity (%)', fontsize=12)
    ax1.set_ylabel('Latent Heat Flux (W/m²)', fontsize=12, color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.axvline(x=base_rh, color='k', linestyle='--', alpha=0.5, label='Monterey Bay RH')
    
    # Mark Monterey Bay condition
    ax1.scatter(base_rh, monterey_latent, color='blue', s=100, marker='*')
    
    # Create a second y-axis for evaporation rate
    ax2 = ax1.twinx()
    ax2.plot(humidity, evaporation, 'r-', linewidth=2)
    ax2.set_ylabel('Evaporation Rate (mm/day)', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Mark Monterey Bay condition
    ax2.scatter(base_rh, monterey_evap, color='red', s=100, marker='*')
    
    # Add title
    plt.title('Effect of Relative Humidity on Latent Heat Flux (Monterey Bay)', fontsize=14)
    plt.grid(False)
    
    # Add a text box with fixed parameters
    textstr = f"Monterey Bay Parameters:\nWind Speed: {wind_speed} m/s\nSea Temp: {sea_temp}°C\nAir Temp: {air_temp}°C"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Add a legend
    ax1.legend(['Latent Heat Flux', 'Monterey Bay RH'], loc='upper right')
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('figures/general_airsea/figures/monterey/monterey_humidity_effect.png', dpi=300)
    plt.close()
    
    print("Humidity effect exploration completed and saved to 'figures/general_airsea/figures/monterey/monterey_humidity_effect.png'")

def compare_airsea_gsw():
    """Compare air-sea flux calculations with GSW density calculations for Monterey Bay water."""
    
    # Create a range of temperatures around Monterey Bay SST
    temp = np.linspace(10, 20, 100)  # 10 to 20°C (range around Monterey Bay)
    
    # Fixed parameters based on Monterey Bay
    salinity = 33.5     # PSU - typical for Monterey Bay
    pressure = 0.0      # dbar (surface)
    latitude = 36.8     # degrees north (Monterey Bay)
    
    # Calculate seawater density using GSW
    # Using SA directly since SA_from_SP requires more parameters than we have
    SA = 33.66    # Approximate Absolute salinity for SP=33.5 at Monterey Bay
    rho_gsw = np.array([gsw.rho(SA, t, pressure) for t in temp])
    
    # Calculate air density for comparison at different temperatures
    humidity = MONTEREY_RH
    air_pressure = MONTEREY_PRESSURE
    rho_air = np.array([atm.air_dens(t, humidity, air_pressure) for t in temp])
    
    # Calculate density ratio (air/water)
    density_ratio = rho_air / rho_gsw
    
    # Mark Monterey Bay condition
    monterey_rho_sw = gsw.rho(SA, MONTEREY_SST, pressure)
    monterey_rho_air = atm.air_dens(MONTEREY_AIR_TEMP, humidity, air_pressure)
    monterey_ratio = monterey_rho_air / monterey_rho_sw
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot water and air density
    ax1.plot(temp, rho_gsw, 'b-', linewidth=2, label='Seawater (GSW)')
    ax1.plot(temp, rho_air, 'r-', linewidth=2, label='Air')
    
    # Mark Monterey Bay conditions
    ax1.scatter(MONTEREY_SST, monterey_rho_sw, color='blue', s=100, marker='*')
    ax1.scatter(MONTEREY_AIR_TEMP, monterey_rho_air, color='red', s=100, marker='*')
    ax1.axvline(x=MONTEREY_SST, color='b', linestyle='--', alpha=0.5, label='Monterey Bay SST')
    ax1.axvline(x=MONTEREY_AIR_TEMP, color='r', linestyle='--', alpha=0.5, label='Monterey Bay Air T')
    
    # Add labels for first subplot
    ax1.set_xlabel('Temperature (°C)', fontsize=12)
    ax1.set_ylabel('Density (kg/m³)', fontsize=12)
    ax1.set_title('Density of Seawater and Air vs. Temperature (Monterey Bay)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Plot density ratio
    ax2.plot(temp, density_ratio * 1000, 'g-', linewidth=2)  # Multiply by 1000 for better visualization
    
    # Mark Monterey Bay condition
    ax2.scatter(MONTEREY_SST, monterey_ratio * 1000, color='green', s=100, marker='*')
    ax2.axvline(x=MONTEREY_SST, color='g', linestyle='--', alpha=0.5, label='Monterey Bay SST')
    
    # Add labels for second subplot
    ax2.set_xlabel('Temperature (°C)', fontsize=12)
    ax2.set_ylabel('Density Ratio (Air/Water) x 1000', fontsize=12)
    ax2.set_title('Density Ratio (Air/Water) vs. Temperature', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend(['Density Ratio', 'Monterey Bay SST'], fontsize=10)
    
    # Add a text box with fixed parameters
    textstr = f"Monterey Bay Parameters:\nSST: {MONTEREY_SST}°C\nAir Temp: {MONTEREY_AIR_TEMP}°C\nSalinity: {salinity} PSU\nRH: {humidity}%"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='bottom', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('figures/general_airsea/figures/monterey/monterey_density_comparison.png', dpi=300)
    plt.close()
    
    print("Density comparison completed and saved to 'figures/general_airsea/figures/monterey/monterey_density_comparison.png'")

def explore_wind_speed_profile():
    """Explore the wind speed profile with height for Monterey Bay."""
    
    # Create a range of heights
    height = np.linspace(1, 100, 100)  # 1 to 100 meters
    
    # Parameters for Monterey Bay
    # Estimate friction velocity from wind speed and roughness length
    u_ref = MONTEREY_WIND_SPEED
    z_ref = 10.0  # Reference height (m)
    z0 = 0.0002    # Roughness length for ocean surface (m)
    
    # Von Karman constant
    kappa = 0.4
    
    # Calculate friction velocity (u*) from log profile
    u_star = u_ref * kappa / np.log(z_ref / z0)
    
    # Calculate wind speed profile using log law
    u = (u_star / kappa) * np.log(height / z0)
    
    # Create the plot
    plt.figure(figsize=(8, 10))
    plt.plot(u, height, 'b-', linewidth=2)
    
    # Add reference heights
    reference_heights = [2, 10, 20, 50]
    for h in reference_heights:
        idx = np.abs(height - h).argmin()
        plt.plot(u[idx], height[idx], 'ro')
        plt.text(u[idx]+0.5, height[idx], f'{h}m: {u[idx]:.1f} m/s', fontsize=10)
    
    # Add labels and title
    plt.xlabel('Wind Speed (m/s)', fontsize=12)
    plt.ylabel('Height (m)', fontsize=12)
    plt.title('Wind Speed Profile with Height over Monterey Bay', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Add a text box explaining the log law and parameters
    textstr = f"Log Law:\nu(z) = (u*/κ) ln(z/z₀)\nwhere:\nu* = {u_star:.2f} m/s (friction velocity)\nκ = 0.4 (von Karman constant)\nz₀ = {z0*1000:.2f} mm (roughness length)"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Add a text box with Monterey Bay parameters
    monterey_str = f"Monterey Bay Parameters:\nWind speed at 10m: {u_ref} m/s"
    props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.5)
    plt.text(0.05, 0.75, monterey_str, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('figures/monterey/monterey_wind_profile.png', dpi=300)
    plt.close()
    
    print("Wind profile exploration completed and saved to 'figures/monterey/monterey_wind_profile.png'")

def main():
    """Run all explorations for Monterey Bay."""
    print("\nStarting air-sea flux exploration for Monterey Bay using the airsea package...")
    print("=" * 75)
    print(f"Using Monterey Bay parameters:")
    print(f"  Wind speed: {MONTEREY_WIND_SPEED:.1f} m/s")
    print(f"  Sea surface temperature: {MONTEREY_SST:.1f} °C")
    print(f"  Air temperature: {MONTEREY_AIR_TEMP:.1f} °C")
    print(f"  Relative humidity: {MONTEREY_RH:.1f} %")
    print(f"  Atmospheric pressure: {MONTEREY_PRESSURE:.1f} mb")
    print("=" * 75)
    
    # Run all exploration functions
    explore_wind_stress_formulations()
    explore_heat_fluxes()
    explore_humidity_effect()
    try:
        compare_airsea_gsw()
    except Exception as e:
        print(f"Error in GSW comparison: {e}")
        print("Skipping GSW comparison...")
    explore_wind_speed_profile()
    
    print("\nAll Monterey Bay explorations complete!")

if __name__ == "__main__":
    main() 
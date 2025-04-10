#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Exploration of air-sea fluxes using the airsea package.
This script examines different flux parameterizations and creates visualizations
to understand relationships between oceanographic and atmospheric variables.
"""

import numpy as np
import matplotlib.pyplot as plt
from airsea import windstress as ws
from airsea import atmosphere as atm
from airsea import constants as const
import gsw  # For comparison with GSW density calculations

# Set plotting style
plt.style.use('ggplot')

def explore_wind_stress_formulations():
    """Compare different wind stress formulations."""
    
    # Create a range of wind speeds to explore (start from 0.1 to avoid divide by zero)
    wind_speeds = np.linspace(0.1, 25, 100)  # 0.1 to 25 m/s
    
    # Calculate wind stress using different formulations
    stress_lp = ws.stress(wind_speeds, drag='largepond')      # Large and Pond (1981)
    stress_s = ws.stress(wind_speeds, drag='smith')           # Smith (1988)
    # Not using 'kondo' as it appears to have issues
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(wind_speeds, stress_lp, 'b-', linewidth=2, label='Large & Pond (1981)')
    plt.plot(wind_speeds, stress_s, 'r-', linewidth=2, label='Smith (1988)')
    
    # Add labels and title
    plt.xlabel('Wind Speed (m/s)', fontsize=12)
    plt.ylabel('Wind Stress (N/m²)', fontsize=12)
    plt.title('Wind Stress vs. Wind Speed: Comparison of Different Formulations', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('wind_stress_comparison.png', dpi=300)
    plt.close()
    
    print("Wind stress comparison completed and saved to 'wind_stress_comparison.png'")

# Custom implementation of sensible heat flux calculation
def calc_sensible_heat(u, Ts, Ta, Pa=1013.0):
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
    Pa : float or array, optional
        Atmospheric pressure (mb), default is 1013.0
        
    Returns:
    --------
    float or array
        Sensible heat flux (W/m²)
    """
    # Air density
    rho_a = atm.air_dens(Ta, 80, Pa)  # Using 80% relative humidity as default
    
    # Specific heat capacity of air
    Cp = 1004.0  # J/(kg·K)
    
    # Stanton number (heat transfer coefficient)
    # Using a simple approximation
    Ch = 1.3e-3
    
    # Calculate sensible heat flux
    H = rho_a * Cp * Ch * u * (Ts - Ta)
    
    return H

# Custom implementation of latent heat flux calculation
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
    
    # Calculate latent heat flux
    E = rho_a * Lv * Ce * u * (q_s - q_a)
    
    return E

def explore_heat_fluxes():
    """Explore sensible and latent heat fluxes under different conditions."""
    
    # Create a range of temperature differences
    delta_t = np.linspace(-5, 5, 100)  # From -5°C to +5°C difference
    
    # Fixed parameters
    wind_speed = 7.0  # m/s
    sea_temp = 15.0   # °C
    air_temp = sea_temp + delta_t  # Varies
    pressure = 1013.0  # hPa
    humidity = 80.0    # Relative humidity (%)
    
    # Calculate heat fluxes
    sensible_heat = np.array([calc_sensible_heat(wind_speed, sea_temp, t) for t in air_temp])
    latent_heat = np.array([calc_latent_heat(wind_speed, sea_temp, t, humidity) for t in air_temp])
    
    # Total heat flux
    total_heat = sensible_heat + latent_heat
    
    # Create the plot
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # Plot sensible heat flux
    ax1.plot(delta_t, sensible_heat, 'b-', linewidth=2, label='Sensible Heat Flux')
    # Plot latent heat flux
    ax1.plot(delta_t, latent_heat, 'r-', linewidth=2, label='Latent Heat Flux')
    # Plot total heat flux
    ax1.plot(delta_t, total_heat, 'g-', linewidth=2, label='Total Heat Flux')
    
    # Add labels and title
    ax1.set_xlabel('Temperature Difference (Sea - Air) °C', fontsize=12)
    ax1.set_ylabel('Heat Flux (W/m²)', fontsize=12)
    ax1.set_title('Heat Fluxes vs. Air-Sea Temperature Difference', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)
    
    # Add a text box explaining sign convention
    textstr = "Positive flux: from ocean to atmosphere (ocean cooling)\nNegative flux: from atmosphere to ocean (ocean warming)"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('heat_flux_comparison.png', dpi=300)
    plt.close()
    
    print("Heat flux comparison completed and saved to 'heat_flux_comparison.png'")

def explore_humidity_effect():
    """Explore the effect of relative humidity on latent heat flux."""
    
    # Create a range of relative humidity values
    humidity = np.linspace(40, 100, 100)  # From 40% to 100%
    
    # Fixed parameters
    wind_speed = 7.0  # m/s
    sea_temp = 20.0   # °C
    air_temp = 18.0   # °C
    pressure = 1013.0  # hPa
    
    # Calculate latent heat flux
    latent_heat = np.array([calc_latent_heat(wind_speed, sea_temp, air_temp, h) for h in humidity])
    
    # Calculate evaporation rate (mm/day)
    evaporation = latent_heat / (2.5e6) * 86400  # Convert W/m² to mm/day using latent heat of vaporization
    
    # Create the plot with two y-axes
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Plot latent heat flux
    ax1.plot(humidity, latent_heat, 'b-', linewidth=2)
    ax1.set_xlabel('Relative Humidity (%)', fontsize=12)
    ax1.set_ylabel('Latent Heat Flux (W/m²)', fontsize=12, color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    
    # Create a second y-axis for evaporation rate
    ax2 = ax1.twinx()
    ax2.plot(humidity, evaporation, 'r-', linewidth=2)
    ax2.set_ylabel('Evaporation Rate (mm/day)', fontsize=12, color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Add title
    plt.title('Effect of Relative Humidity on Latent Heat Flux and Evaporation', fontsize=14)
    plt.grid(False)
    
    # Add a text box with fixed parameters
    textstr = f"Fixed Parameters:\nWind Speed: {wind_speed} m/s\nSea Temp: {sea_temp}°C\nAir Temp: {air_temp}°C"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('humidity_effect.png', dpi=300)
    plt.close()
    
    print("Humidity effect exploration completed and saved to 'humidity_effect.png'")

def compare_airsea_gsw():
    """Compare air-sea flux calculations with GSW density calculations."""
    
    # Create a range of temperatures
    temp = np.linspace(0, 30, 100)  # 0 to 30°C
    
    # Fixed parameters
    salinity = 35.0     # PSU
    pressure = 0.0      # dbar (surface)
    latitude = 45.0     # degrees north
    
    # Calculate seawater density using GSW
    # Using SA directly since SA_from_SP requires more parameters than we have
    SA = 35.16504    # Absolute salinity for SP=35 at lat=45, lon=0, p=0
    rho_gsw = np.array([gsw.rho(SA, t, pressure) for t in temp])
    
    # Calculate air density for comparison at different temperatures
    humidity = 80.0    # %
    air_pressure = 1013.0  # hPa
    rho_air = np.array([atm.air_dens(t, humidity, air_pressure) for t in temp])
    
    # Calculate density ratio (air/water)
    density_ratio = rho_air / rho_gsw
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot water and air density
    ax1.plot(temp, rho_gsw, 'b-', linewidth=2, label='Seawater (GSW)')
    ax1.plot(temp, rho_air, 'r-', linewidth=2, label='Air')
    
    # Add labels for first subplot
    ax1.set_xlabel('Temperature (°C)', fontsize=12)
    ax1.set_ylabel('Density (kg/m³)', fontsize=12)
    ax1.set_title('Density of Seawater and Air vs. Temperature', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)
    
    # Plot density ratio
    ax2.plot(temp, density_ratio * 1000, 'g-', linewidth=2)  # Multiply by 1000 for better visualization
    
    # Add labels for second subplot
    ax2.set_xlabel('Temperature (°C)', fontsize=12)
    ax2.set_ylabel('Density Ratio (Air/Water) x 1000', fontsize=12)
    ax2.set_title('Density Ratio (Air/Water) vs. Temperature', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # Add a text box with fixed parameters
    textstr = f"Fixed Parameters:\nSalinity: {salinity} PSU\nAir Pressure: {air_pressure} hPa\nRelative Humidity: {humidity}%"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.05, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='bottom', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('density_comparison.png', dpi=300)
    plt.close()
    
    print("Density comparison completed and saved to 'density_comparison.png'")

def custom_drag_coefficient(wind_speed, formulation='largepond'):
    """
    Calculate the drag coefficient based on different formulations.
    
    Parameters:
    -----------
    wind_speed : float or array
        Wind speed (m/s)
    formulation : str, optional
        Formulation to use ('largepond' or 'smith'), default is 'largepond'
        
    Returns:
    --------
    float or array
        Drag coefficient
    """
    if formulation == 'largepond':
        # Large and Pond (1981) formulation
        cd = np.ones_like(wind_speed) * 1.2e-3  # Initial value
        
        # Adjust for wind speeds > 11 m/s
        mask = wind_speed > 11
        cd[mask] = 0.49e-3 + 0.065e-3 * wind_speed[mask]
        
        return cd
        
    elif formulation == 'smith':
        # Smith (1988) formulation
        a = 0.61 + 0.063 * wind_speed
        cd = np.ones_like(wind_speed) * 1.0e-3  # Initial value
        cd = (a / 1000) * np.ones_like(wind_speed)
        
        return cd
    else:
        raise ValueError(f"Unknown formulation: {formulation}")

def explore_drag_coefficient():
    """Explore how drag coefficient varies with wind speed for different formulations."""
    
    # Create a range of wind speeds (start from 0.1 to avoid divide by zero)
    wind_speeds = np.linspace(0.1, 30, 100)  # 0.1 to 30 m/s
    
    # Calculate drag coefficients using our custom function
    Cd_lp = custom_drag_coefficient(wind_speeds, formulation='largepond')
    Cd_smith = custom_drag_coefficient(wind_speeds, formulation='smith')
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(wind_speeds, Cd_lp * 1000, 'b-', linewidth=2, label='Large & Pond (1981)')
    plt.plot(wind_speeds, Cd_smith * 1000, 'r-', linewidth=2, label='Smith (1988)')
    
    # Add labels and title
    plt.xlabel('Wind Speed (m/s)', fontsize=12)
    plt.ylabel('Drag Coefficient (Cd × 1000)', fontsize=12)
    plt.title('Drag Coefficient vs. Wind Speed: Comparison of Different Formulations', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('drag_coefficient_comparison.png', dpi=300)
    plt.close()
    
    print("Drag coefficient comparison completed and saved to 'drag_coefficient_comparison.png'")

def explore_wind_speed_profile():
    """Explore the wind speed profile with height."""
    
    # Create a range of heights
    height = np.linspace(1, 100, 100)  # 1 to 100 meters
    
    # Fixed parameters
    u_star = 0.5  # friction velocity (m/s)
    z0 = 0.001    # roughness length (m)
    
    # Calculate wind speed profile using log law
    kappa = 0.4  # von Karman constant
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
    plt.title('Wind Speed Profile with Height (Logarithmic Law)', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Add a text box explaining the log law
    textstr = "Log Law:\nu(z) = (u*/κ) ln(z/z₀)\nwhere:\nu* = friction velocity\nκ = von Karman constant (0.4)\nz₀ = roughness length"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top', bbox=props)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig('wind_profile.png', dpi=300)
    plt.close()
    
    print("Wind profile exploration completed and saved to 'wind_profile.png'")

def main():
    """Run all explorations."""
    print("\nStarting air-sea flux exploration using the airsea package...")
    
    # Run all exploration functions
    explore_wind_stress_formulations()
    explore_heat_fluxes()
    explore_humidity_effect()
    try:
        compare_airsea_gsw()
    except Exception as e:
        print(f"Error in GSW comparison: {e}")
        print("Skipping GSW comparison...")
    explore_drag_coefficient()
    explore_wind_speed_profile()
    
    print("\nAll explorations complete!")

if __name__ == "__main__":
    main() 
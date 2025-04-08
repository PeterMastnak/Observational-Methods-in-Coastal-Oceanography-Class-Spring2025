#!/usr/bin/env python
"""
Example script demonstrating the use of the GSW package for oceanographic calculations.
This script shows basic computations such as density, buoyancy frequency, and mixing parameters.
"""

import numpy as np
import matplotlib.pyplot as plt
import gsw
import cmocean  # Specialized colormaps for oceanography

def main():
    print("GSW-Python Examples for Coastal Oceanography")
    print("--------------------------------------------")
    
    # Example 1: Calculate seawater density from practical salinity, temperature, and pressure
    # These could be data from a CTD cast (Conductivity, Temperature, Depth)
    print("\nExample 1: Calculate seawater density")
    
    # Sample data (typical values from coastal regions)
    # Practical salinity (PSU) - unitless 
    SP = np.array([35.5, 35.7, 36.0, 36.2, 36.4])
    # Temperature (°C)
    t = np.array([20.0, 18.0, 15.0, 10.0, 8.0])
    # Pressure (dbar) - approximately equal to depth in meters
    p = np.array([0, 200, 400, 600, 800])
    
    # Location for calculations
    lon = -70  # Longitude
    lat = 40   # Latitude
    
    # Convert practical salinity to absolute salinity
    # Note: In gsw, lon and lat are positional arguments, not keyword arguments
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    print(f"Practical Salinity (PSU): {SP}")
    print(f"Absolute Salinity (g/kg): {SA.round(3)}")
    
    # Convert temperature to conservative temperature
    CT = gsw.CT_from_t(SA, t, p)
    print(f"In-situ Temperature (°C): {t}")
    print(f"Conservative Temperature (°C): {CT.round(3)}")
    
    # Calculate density
    rho = gsw.rho(SA, CT, p)
    print(f"Seawater Density (kg/m³): {rho.round(3)}")
    
    # Potential density anomaly (sigma-0)
    sigma0 = gsw.sigma0(SA, CT)
    print(f"Potential Density Anomaly (kg/m³): {sigma0.round(3)}")
    
    # Example 2: Calculate the buoyancy (Brunt-Väisälä) frequency
    print("\nExample 2: Calculate buoyancy frequency")
    N2, p_mid = gsw.Nsquared(SA, CT, p, lat)
    print(f"Buoyancy Frequency Squared (rad²/s²): {N2.round(6)}")
    print(f"At pressure levels (dbar): {p_mid}")
    
    # Calculate buoyancy frequency in cycles per hour
    N_cph = np.sqrt(N2) * 3600 / (2 * np.pi)
    print(f"Buoyancy Frequency (cycles per hour): {N_cph.round(2)}")
    
    # Example 3: T-S Diagram 
    print("\nExample 3: Creating a T-S Diagram")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create sigma-0 contours for the T-S diagram background
    SA_grid = np.linspace(34.0, 37.0, 100)
    CT_grid = np.linspace(5.0, 25.0, 100)
    SA_grid, CT_grid = np.meshgrid(SA_grid, CT_grid)
    sigma0_grid = gsw.sigma0(SA_grid, CT_grid)
    
    # Plot density contours
    CS = ax.contour(SA_grid, CT_grid, sigma0_grid, colors='grey', alpha=0.5)
    ax.clabel(CS, inline=1, fontsize=10, fmt='%1.1f')
    
    # Plot the actual data
    sc = ax.scatter(SA, CT, c=p, cmap=cmocean.cm.deep, s=80, label='Sample Data')
    ax.plot(SA, CT, 'k--', alpha=0.5)
    
    # Labels and title
    ax.set_xlabel('Absolute Salinity (g/kg)')
    ax.set_ylabel('Conservative Temperature (°C)')
    ax.set_title('Temperature-Salinity (T-S) Diagram with σ₀ Contours')
    
    # Add colorbar for pressure
    cbar = plt.colorbar(sc)
    cbar.set_label('Pressure (dbar)')
    
    # Save the figure
    plt.savefig('ts_diagram.png', dpi=300, bbox_inches='tight')
    print("T-S diagram saved as 'ts_diagram.png'")
    
    print("\nExample 4: Water mass mixing")
    # Define two water masses
    SA1, CT1 = 35.0, 18.0  # Water mass 1
    SA2, CT2 = 36.5, 10.0  # Water mass 2
    
    # Calculate mixing proportions
    mixing_proportions = np.linspace(0, 1, 11)
    SA_mix = SA1 + mixing_proportions * (SA2 - SA1)
    CT_mix = CT1 + mixing_proportions * (CT2 - CT1)
    
    # Calculate density of the mixed water
    rho_mix = gsw.rho(SA_mix, CT_mix, 0)
    
    # Calculate density if simple mixing (linear)
    rho1 = gsw.rho(SA1, CT1, 0)
    rho2 = gsw.rho(SA2, CT2, 0)
    rho_linear = rho1 + mixing_proportions * (rho2 - rho1)
    
    # Calculate the difference (cabbeling effect)
    cabbeling = rho_mix - rho_linear
    
    print(f"Mixing Proportions: {mixing_proportions}")
    print(f"Mixed Water Density (kg/m³): {rho_mix.round(3)}")
    print(f"Linear Mixing Density (kg/m³): {rho_linear.round(3)}")
    print(f"Cabbeling Effect (kg/m³): {cabbeling.round(3)}")
    
    print("\nComplete! Run this script with matplotlib interactive mode to view the plots.")

if __name__ == "__main__":
    main() 
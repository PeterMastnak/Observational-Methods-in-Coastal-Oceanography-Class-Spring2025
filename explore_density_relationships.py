#!/usr/bin/env python
"""
Explore the effects of pressure, salinity, and temperature on seawater density.
This script also compares the differences between Absolute Salinity vs. Practical Salinity
and Conservative Temperature vs. In situ Temperature.

This is a Python implementation of the group exercise that uses the GSW package
(Python version of the TEOS-10 Gibbs SeaWater Oceanographic Toolbox).
"""

import numpy as np
import matplotlib.pyplot as plt
import gsw
from matplotlib.gridspec import GridSpec
import cmocean

def explore_density_pressure_relationship():
    """Explore how pressure affects density with fixed temperature and salinity."""
    # Fixed parameters
    SP = 35.0  # Practical Salinity (PSU)
    t = 10.0   # Temperature (°C)
    lon, lat = -70, 40  # Approximate location (Western Atlantic)
    
    # Typical ocean depth range (0 to 10000 meters, roughly equivalent to 0-10000 dbar)
    p_range = np.linspace(0, 10000, 100)  # Pressure in dbar
    
    # Convert practical salinity to absolute salinity
    SA = gsw.SA_from_SP(SP, p_range, lon, lat)
    
    # Convert in-situ temperature to conservative temperature
    CT = gsw.CT_from_t(SA, t, p_range)
    
    # Calculate density at each pressure
    rho = gsw.rho(SA, CT, p_range)
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(p_range, rho, 'b-', linewidth=2)
    ax.set_xlabel('Pressure (dbar)', fontsize=12)
    ax.set_ylabel('Density (kg/m³)', fontsize=12)
    ax.set_title(f'Effect of Pressure on Density\n(Fixed Temperature: {t}°C, Practical Salinity: {SP})', fontsize=14)
    ax.grid(True)
    
    # Calculate and display the density change
    density_change = rho[-1] - rho[0]
    percent_change = (density_change / rho[0]) * 100
    ax.text(0.05, 0.95, f"Density at surface: {rho[0]:.2f} kg/m³\nDensity at 10000m: {rho[-1]:.2f} kg/m³\n"
             f"Change: {density_change:.2f} kg/m³ ({percent_change:.2f}%)",
             transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=12)
    
    plt.tight_layout()
    plt.savefig('density_vs_pressure.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'density_vs_pressure.png'")
    
    return fig, ax

def explore_density_temperature_relationship():
    """Explore how temperature affects density with fixed pressure and salinity."""
    # Fixed parameters
    SP = 35.0  # Practical Salinity (PSU)
    p = 0.0    # Surface pressure (dbar)
    lon, lat = -70, 40  # Approximate location (Western Atlantic)
    
    # Typical ocean temperature range
    t_range = np.linspace(-2, 30, 100)  # Temperature in °C
    
    # Convert practical salinity to absolute salinity (constant at fixed pressure)
    SA = gsw.SA_from_SP(SP, p, lon, lat) * np.ones_like(t_range)
    
    # Convert in-situ temperature to conservative temperature
    CT = gsw.CT_from_t(SA, t_range, p)
    
    # Calculate density at each temperature
    rho = gsw.rho(SA, CT, p)
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_range, rho, 'r-', linewidth=2)
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Density (kg/m³)', fontsize=12)
    ax.set_title(f'Effect of Temperature on Density\n(Fixed Pressure: {p} dbar, Practical Salinity: {SP})', fontsize=14)
    ax.grid(True)
    
    # Calculate and display the density change
    density_change = np.max(rho) - np.min(rho)
    percent_change = (density_change / np.max(rho)) * 100
    ax.text(0.05, 0.05, f"Density range: {np.min(rho):.2f} to {np.max(rho):.2f} kg/m³\n"
             f"Change: {density_change:.2f} kg/m³ ({percent_change:.2f}%)",
             transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=12)

    # Find temperature of maximum density (if it exists in our range)
    max_density_idx = np.argmax(rho)
    max_density_temp = t_range[max_density_idx]
    if max_density_idx > 0 and max_density_idx < len(t_range) - 1:
        ax.axvline(x=max_density_temp, color='k', linestyle='--', alpha=0.5)
        ax.text(max_density_temp + 1, np.max(rho), f'Tₘₐₓ ≈ {max_density_temp:.2f}°C',
                va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('density_vs_temperature.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'density_vs_temperature.png'")
    
    return fig, ax

def explore_density_salinity_relationship():
    """Explore how salinity affects density with fixed pressure and temperature."""
    # Fixed parameters
    t = 10.0    # Temperature (°C)
    p = 0.0     # Surface pressure (dbar)
    lon, lat = -70, 40  # Approximate location (Western Atlantic)
    
    # Typical ocean salinity range
    SP_range = np.linspace(0, 40, 100)  # Practical Salinity (PSU)
    
    # Convert practical salinity to absolute salinity
    SA_range = np.array([gsw.SA_from_SP(sp, p, lon, lat) for sp in SP_range])
    
    # Convert in-situ temperature to conservative temperature (constant)
    CT = gsw.CT_from_t(SA_range, t, p)
    
    # Calculate density at each salinity
    rho = gsw.rho(SA_range, CT, p)
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(SP_range, rho, 'g-', linewidth=2)
    ax.set_xlabel('Practical Salinity (PSU)', fontsize=12)
    ax.set_ylabel('Density (kg/m³)', fontsize=12)
    ax.set_title(f'Effect of Salinity on Density\n(Fixed Temperature: {t}°C, Pressure: {p} dbar)', fontsize=14)
    ax.grid(True)
    
    # Calculate and display the density change
    density_change = rho[-1] - rho[0]
    percent_change = (density_change / rho[0]) * 100
    ax.text(0.05, 0.05, f"Fresh water density: {rho[0]:.2f} kg/m³\nSalty water density: {rho[-1]:.2f} kg/m³\n"
             f"Change: {density_change:.2f} kg/m³ ({percent_change:.2f}%)",
             transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=12)
    
    # Mark typical ocean salinity of ~35 PSU
    typical_idx = np.abs(SP_range - 35).argmin()
    ax.axvline(x=SP_range[typical_idx], color='k', linestyle='--', alpha=0.5)
    ax.text(SP_range[typical_idx] + 0.5, rho[typical_idx], f'Typical Ocean Salinity ≈ 35 PSU',
            va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('density_vs_salinity.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'density_vs_salinity.png'")
    
    return fig, ax

def compare_salinity_measures():
    """Compare Absolute Salinity vs Practical Salinity."""
    # Fixed parameters
    p = 0.0  # Surface pressure (dbar)
    
    # Create a grid of longitudes and latitudes
    lon = np.linspace(-180, 180, 100)
    lat = np.linspace(-80, 80, 80)
    LON, LAT = np.meshgrid(lon, lat)
    
    # Use a constant Practical Salinity
    SP = 35.0 * np.ones_like(LON)
    P = p * np.ones_like(LON)
    
    # Calculate Absolute Salinity
    SA = gsw.SA_from_SP(SP, P, LON, LAT)
    
    # Calculate the difference
    SA_diff = SA - SP
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(12, 8))
    contour = ax.contourf(LON, LAT, SA_diff, cmap=cmocean.cm.haline, levels=20)
    ax.set_xlabel('Longitude (°E)', fontsize=12)
    ax.set_ylabel('Latitude (°N)', fontsize=12)
    ax.set_title('Difference Between Absolute and Practical Salinity (g/kg - PSU)\n'
                 'at Fixed Practical Salinity of 35 PSU', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('SA - SP (g/kg - PSU)', fontsize=12)
    
    # Add annotations
    max_diff = np.max(SA_diff)
    min_diff = np.min(SA_diff)
    max_loc = np.unravel_index(np.argmax(SA_diff), SA_diff.shape)
    min_loc = np.unravel_index(np.argmin(SA_diff), SA_diff.shape)
    
    ax.text(0.02, 0.02, f"Range: {min_diff:.3f} to {max_diff:.3f} g/kg - PSU\n"
                       f"Mean difference: {np.mean(SA_diff):.3f} g/kg - PSU",
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=10)
    
    plt.tight_layout()
    plt.savefig('absolute_vs_practical_salinity.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'absolute_vs_practical_salinity.png'")
    
    return fig, ax

def compare_temperature_measures():
    """Compare Conservative Temperature vs In-situ Temperature."""
    # Fixed parameters
    SP = 35.0  # Practical Salinity (PSU)
    lon, lat = -70, 40  # Approximate location (Western Atlantic)
    
    # Create a grid of temperatures and pressures
    t = np.linspace(-2, 30, 50)  # In-situ temperature in °C
    p = np.linspace(0, 5000, 50)  # Pressure in dbar
    T, P = np.meshgrid(t, p)
    
    # Convert practical salinity to absolute salinity
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    
    # Convert in-situ temperature to conservative temperature
    CT = gsw.CT_from_t(SA, T, P)
    
    # Calculate the difference
    temp_diff = CT - T
    
    # Create a figure
    fig, ax = plt.subplots(figsize=(12, 8))
    contour = ax.contourf(T, P, temp_diff, cmap=cmocean.cm.thermal, levels=20)
    ax.set_xlabel('In-situ Temperature (°C)', fontsize=12)
    ax.set_ylabel('Pressure (dbar)', fontsize=12)
    ax.set_title('Difference Between Conservative and In-situ Temperature (°C)\n'
                 f'at Fixed Practical Salinity of {SP} PSU', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Invert y-axis to show depth increasing downward
    ax.invert_yaxis()
    
    # Add colorbar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('CT - t (°C)', fontsize=12)
    
    # Add annotations
    max_diff = np.max(temp_diff)
    min_diff = np.min(temp_diff)
    
    ax.text(0.02, 0.02, f"Range: {min_diff:.3f} to {max_diff:.3f} °C\n"
                       f"Maximum difference at high temperature and pressure",
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7), fontsize=10)
    
    plt.tight_layout()
    plt.savefig('conservative_vs_insitu_temperature.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'conservative_vs_insitu_temperature.png'")
    
    return fig, ax

def density_3d_visualization():
    """Create a 3D visualization of how temperature, salinity, and pressure affect density."""
    # Create a grid of temperatures and salinities at two different pressures
    t = np.linspace(0, 30, 30)  # Temperature in °C
    SP = np.linspace(30, 40, 30)  # Practical Salinity
    T, SP_grid = np.meshgrid(t, SP)
    
    # Location for calculating Absolute Salinity
    lon, lat = -70, 40
    
    # Calculate density at surface (0 dbar)
    p_surface = 0
    SA_surface = gsw.SA_from_SP(SP_grid, p_surface, lon, lat)
    CT_surface = gsw.CT_from_t(SA_surface, T, p_surface)
    rho_surface = gsw.rho(SA_surface, CT_surface, p_surface)
    
    # Calculate density at 1000m depth (approximately 1000 dbar)
    p_deep = 1000
    SA_deep = gsw.SA_from_SP(SP_grid, p_deep, lon, lat)
    CT_deep = gsw.CT_from_t(SA_deep, T, p_deep)
    rho_deep = gsw.rho(SA_deep, CT_deep, p_deep)
    
    # Create a figure with two subplots using GridSpec
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], figure=fig)
    
    # Surface density plot (0 dbar)
    ax1 = fig.add_subplot(gs[0, 0])
    cs1 = ax1.contourf(T, SP_grid, rho_surface, levels=20, cmap=cmocean.cm.dense)
    ax1.set_xlabel('Temperature (°C)', fontsize=12)
    ax1.set_ylabel('Practical Salinity (PSU)', fontsize=12)
    ax1.set_title('Density (kg/m³) at Surface (0 dbar)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Deep density plot (1000 dbar)
    ax2 = fig.add_subplot(gs[0, 1])
    cs2 = ax2.contourf(T, SP_grid, rho_deep, levels=20, cmap=cmocean.cm.dense)
    ax2.set_xlabel('Temperature (°C)', fontsize=12)
    ax2.set_ylabel('Practical Salinity (PSU)', fontsize=12)
    ax2.set_title('Density (kg/m³) at 1000 dbar', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # Add colorbar
    cax = fig.add_subplot(gs[0, 2])
    cbar = plt.colorbar(cs2, cax=cax)
    cbar.set_label('Density (kg/m³)', fontsize=12)
    
    # Add density difference statistics
    density_diff = rho_deep - rho_surface
    fig.text(0.5, 0.02, 
             f"Mean density at surface: {np.mean(rho_surface):.2f} kg/m³\n"
             f"Mean density at 1000 dbar: {np.mean(rho_deep):.2f} kg/m³\n"
             f"Mean density increase with depth: {np.mean(density_diff):.2f} kg/m³ ({np.mean(density_diff)/np.mean(rho_surface)*100:.2f}%)",
             ha='center', bbox=dict(facecolor='white', alpha=0.7), fontsize=12)
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig('figures/general_airsea/density_3d_comparison.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'figures/general_airsea/density_3d_comparison.png'")
    
    return fig

def main():
    """Run all the analyses and create summary findings."""
    print("Exploring how pressure, salinity, and temperature affect seawater density...\n")
    
    # Run each analysis
    explore_density_pressure_relationship()
    explore_density_temperature_relationship()
    explore_density_salinity_relationship()
    compare_salinity_measures()
    compare_temperature_measures()
    density_3d_visualization()
    
    # Generate a summary report
    with open('density_relationships_summary.txt', 'w') as f:
        f.write("SUMMARY OF DENSITY RELATIONSHIPS IN SEAWATER\n")
        f.write("===========================================\n\n")
        
        f.write("1. EFFECTS OF PRESSURE ON DENSITY:\n")
        f.write("   - Seawater density increases with pressure (depth).\n")
        f.write("   - For typical ocean parameters, density increases by approximately 4-5% \n")
        f.write("     from the surface to the deep ocean (~10,000m).\n")
        f.write("   - This effect is nearly linear and is due to the compression of water molecules.\n\n")
        
        f.write("2. EFFECTS OF TEMPERATURE ON DENSITY:\n")
        f.write("   - Seawater density generally decreases as temperature increases.\n")
        f.write("   - The relationship is non-linear, with the greatest effect in colder waters.\n")
        f.write("   - Unlike freshwater, seawater's density typically does not have a maximum \n")
        f.write("     at 4°C due to the presence of salt.\n")
        f.write("   - Temperature can cause density variations of approximately 1-2% \n")
        f.write("     across oceanic temperature ranges (typically -2°C to 30°C).\n\n")
        
        f.write("3. EFFECTS OF SALINITY ON DENSITY:\n")
        f.write("   - Seawater density increases with salinity.\n")
        f.write("   - The relationship is nearly linear.\n")
        f.write("   - Across typical ocean salinity ranges (33-37 PSU), density varies by around 0.3-0.4%.\n")
        f.write("   - When comparing fresh water to typical seawater, the difference is about 2-3%.\n\n")
        
        f.write("4. ABSOLUTE VS. PRACTICAL SALINITY:\n")
        f.write("   - Absolute Salinity (SA, g/kg) is the mass fraction of dissolved material in seawater.\n")
        f.write("   - Practical Salinity (SP, PSU) is a conductivity-based scale without units.\n")
        f.write("   - The difference between them varies geographically due to variations in seawater composition.\n")
        f.write("   - The difference is typically 0.1-0.5 g/kg higher for Absolute Salinity.\n")
        f.write("   - Largest differences occur in the North Pacific and Indian Ocean due to silicate content.\n\n")
        
        f.write("5. CONSERVATIVE VS. IN-SITU TEMPERATURE:\n")
        f.write("   - Conservative Temperature (CT) represents the heat content of seawater.\n")
        f.write("   - In-situ Temperature (t) is the measured temperature at a specific pressure.\n")
        f.write("   - CT is more appropriate for heat budget calculations.\n")
        f.write("   - Differences increase with pressure and are largest in warm, deep waters.\n")
        f.write("   - CT is typically slightly lower than in-situ temperature at depth.\n\n")
        
        f.write("6. COMBINED EFFECTS:\n")
        f.write("   - Of the three factors, temperature variations typically cause the largest density changes \n")
        f.write("     in the upper ocean, while pressure dominates in the deep ocean.\n")
        f.write("   - In certain regions like the polar oceans, salinity can be the dominant factor controlling density.\n")
        f.write("   - These density differences drive global ocean circulation patterns.\n")
        f.write("   - The nonlinear nature of the seawater equation of state leads to important phenomena \n")
        f.write("     like cabbeling, where mixing of water masses can increase density.\n\n")
        
        f.write("Generated using GSW-Python (Gibbs SeaWater) toolbox based on TEOS-10.")
    
    print("\nAll analyses complete! Summary saved to 'density_relationships_summary.txt'")
    print("Open the PNG image files to view the visualizations.")

if __name__ == "__main__":
    main() 
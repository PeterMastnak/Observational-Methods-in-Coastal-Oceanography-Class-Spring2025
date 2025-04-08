#!/usr/bin/env python
"""
Example script demonstrating how to read and process oceanographic data from NetCDF files.
This script shows how to access and analyze data typically used in coastal oceanography studies.
"""

import numpy as np
import matplotlib.pyplot as plt
import gsw
import xarray as xr
import cmocean
import os
import urllib.request
from pathlib import Path

def download_sample_data(url, filename):
    """Download sample oceanographic data if it doesn't exist."""
    if not os.path.exists(filename):
        print(f"Downloading sample data from {url}...")
        Path(os.path.dirname(filename)).mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(url, filename)
        print(f"Downloaded to {filename}")
    else:
        print(f"Using existing data file: {filename}")

def main():
    print("Processing Ocean Data Example")
    print("----------------------------")
    
    # Create a data directory if it doesn't exist
    data_dir = 'data'
    os.makedirs(data_dir, exist_ok=True)
    
    # Define the path to our local synthetic data file
    filename = os.path.join(data_dir, 'synthetic_ocean_temp.nc')
    
    # Check if the synthetic data file exists
    if not os.path.exists(filename):
        print(f"Synthetic data file not found at {filename}")
        print("Please run 'python generate_sample_data.py' first to create the synthetic data file.")
        return
    
    try:
        # Read the NetCDF data using xarray
        print("\nReading NetCDF data...")
        print(f"Using synthetic ocean temperature data: {filename}")
        ds = xr.open_dataset(filename)
        
        print("\nDataset overview:")
        print(ds)
        
        # Extract data for analysis
        print("\nExtracting data for analysis...")
        # Get temperature data (first time point)
        temp = ds['t_an'].isel(time=0)
        # Get depth, latitude, and longitude
        depth = ds['depth']
        lat = ds['lat']
        lon = ds['lon']
        
        print(f"Data dimensions: {temp.shape}")
        print(f"Depth levels: {len(depth)} (from {depth.min().item():.1f} to {depth.max().item():.1f} meters)")
        
        # Create a transect along a specific latitude
        target_lat = 40.0  # North Atlantic, approximately
        target_lon_range = (-75, -65)  # US East Coast to mid-Atlantic
        
        # Find closest latitude index
        lat_idx = np.abs(lat.values - target_lat).argmin()
        actual_lat = lat.values[lat_idx]
        print(f"\nCreating transect at latitude {actual_lat:.2f}°N")
        
        # Get longitude indices
        lon_indices = np.where((lon.values >= target_lon_range[0]) & (lon.values <= target_lon_range[1]))[0]
        
        # Check if we have any longitude indices
        if len(lon_indices) == 0:
            # If no exact match, take a range centered on the middle of our longitude array
            mid_lon_idx = len(lon) // 2
            half_width = min(5, len(lon) // 4)  # Use at most 1/4 of the available longitudes
            lon_indices = np.arange(mid_lon_idx - half_width, mid_lon_idx + half_width)
            print(f"No longitudes in the specified range. Using longitudes {lon.values[lon_indices[0]]:.2f}° to {lon.values[lon_indices[-1]]:.2f}°E instead.")
        
        transect_lons = lon.values[lon_indices]
        
        # Extract temperature transect
        temp_transect = temp.isel(lat=lat_idx).isel(lon=lon_indices)
        
        # Create a figure with a cross-section
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Create a mesh for contour plot
        lon_mesh, depth_mesh = np.meshgrid(transect_lons, depth.values)
        
        # Plot temperature cross-section
        cf = ax.contourf(lon_mesh, depth_mesh, temp_transect.values, 
                        levels=20, cmap=cmocean.cm.thermal)
        
        # Add contour lines
        cs = ax.contour(lon_mesh, depth_mesh, temp_transect.values, 
                        levels=np.arange(0, 30, 2), colors='k', linewidths=0.5)
        ax.clabel(cs, inline=1, fontsize=8, fmt='%d')
        
        # Set axis properties
        ax.set_ylim(depth.max(), 0)  # Invert y-axis
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Depth (m)')
        ax.set_title(f'Temperature (°C) at {actual_lat:.2f}°N - Synthetic Ocean Data')
        
        # Add colorbar
        cbar = plt.colorbar(cf, ax=ax)
        cbar.set_label('Temperature (°C)')
        
        # Save the figure
        plt.tight_layout()
        output_file = os.path.join(data_dir, 'temperature_transect.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_file}")
        
        # Create a map view of temperature at specific depth
        # Choose a depth level (e.g., 100m)
        target_depth = 100
        depth_idx = np.abs(depth.values - target_depth).argmin()
        actual_depth = depth.values[depth_idx]
        
        # Extract temperature at this depth
        temp_map = temp.isel(depth=depth_idx)
        
        print(f"\nCreating map view at depth {actual_depth:.1f}m")
        
        # Plot map
        try:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            
            # Plot temperature map
            map_extent = [lon.min().item(), lon.max().item(), 
                         lat.min().item(), lat.max().item()]
            
            fig = plt.figure(figsize=(12, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
            ax.set_extent(map_extent)
            
            # Add map features
            ax.add_feature(cfeature.LAND, facecolor='lightgray')
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS, linestyle=':')
            ax.gridlines(draw_labels=True)
            
            # Plot temperature
            cf = ax.contourf(lon.values, lat.values, temp_map.values, 
                            levels=20, cmap=cmocean.cm.thermal, transform=ccrs.PlateCarree())
            
            # Mark the transect
            ax.plot([transect_lons[0], transect_lons[-1]], [actual_lat, actual_lat], 
                    'r-', linewidth=2, transform=ccrs.PlateCarree())
            
            # Add colorbar
            cbar = plt.colorbar(cf, ax=ax, shrink=0.8)
            cbar.set_label('Temperature (°C)')
            
            # Add title
            ax.set_title(f'Temperature at {actual_depth:.1f}m - Synthetic Ocean Data')
            
            # Save the figure
            output_file = os.path.join(data_dir, 'temperature_map.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Map figure saved to {output_file}")
            
        except ImportError:
            print("Cartopy not properly installed. Skipping map creation.")
            # Create a simple contour plot instead
            fig, ax = plt.subplots(figsize=(10, 8))
            cf = ax.contourf(lon.values, lat.values, temp_map.values, levels=20, cmap=cmocean.cm.thermal)
            ax.set_xlabel('Longitude (°E)')
            ax.set_ylabel('Latitude (°N)')
            ax.set_title(f'Temperature at {actual_depth:.1f}m - Synthetic Ocean Data')
            ax.grid(True)
            
            # Add colorbar
            cbar = plt.colorbar(cf, ax=ax)
            cbar.set_label('Temperature (°C)')
            
            # Mark the transect line
            ax.plot([transect_lons[0], transect_lons[-1]], [actual_lat, actual_lat], 'r-', linewidth=2)
            
            # Save the figure
            output_file = os.path.join(data_dir, 'temperature_map_simple.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Simple map figure saved to {output_file}")
            
        # Close the dataset
        ds.close()
        
        # Use GSW to calculate derived properties from our synthetic data
        print("\nExample calculation with GSW (using our synthetic data):")
        
        # Create artificial salinity data for demonstration
        # A simple salinity profile that increases with depth
        SP = 35.0 + np.linspace(0, 1.5, len(depth))  # Increasing with depth
        p = depth.values  # Pressure approximately equals depth in dbar
        
        # Extract temperature along the center of our domain
        center_lat_idx = len(lat) // 2
        center_lon_idx = len(lon) // 2
        t = temp.isel(lat=center_lat_idx, lon=center_lon_idx).values
        
        # Location for demonstration
        lon_demo, lat_demo = lon.values[center_lon_idx], lat.values[center_lat_idx]
        
        # Convert to absolute salinity and conservative temperature
        SA = gsw.SA_from_SP(SP, p, lon_demo, lat_demo)
        CT = gsw.CT_from_t(SA, t, p)
        
        # Calculate density
        rho = gsw.rho(SA, CT, p)
        
        # Calculate buoyancy frequency (N²)
        N2, p_mid = gsw.Nsquared(SA, CT, p, lat_demo)
        
        # Plot temperature, salinity, density, and N² profiles
        fig, axes = plt.subplots(1, 4, figsize=(15, 8), sharey=True)
        
        # Temperature
        axes[0].plot(t, p, 'r-')
        axes[0].set_xlabel('Temperature (°C)')
        axes[0].set_ylabel('Depth (m)')
        axes[0].set_ylim(max(p), 0)  # Invert y-axis
        axes[0].grid(True)
        
        # Salinity
        axes[1].plot(SP, p, 'b-')
        axes[1].set_xlabel('Practical Salinity')
        axes[1].grid(True)
        
        # Density
        axes[2].plot(rho - 1000, p, 'g-')
        axes[2].set_xlabel('Density Anomaly (kg/m³)')
        axes[2].grid(True)
        
        # N²
        N_cph = np.sqrt(N2) * 3600 / (2 * np.pi)  # Convert to cycles per hour
        axes[3].plot(N_cph, p_mid, 'k-')
        axes[3].set_xlabel('Buoyancy Freq. (cph)')
        axes[3].grid(True)
        
        # Add title
        plt.suptitle('Synthetic Ocean Profiles with GSW Calculations', fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        output_file = os.path.join(data_dir, 'ocean_profiles.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Profile figure saved to {output_file}")
        
        print("\nComplete! Run this script with matplotlib interactive mode to view the plots.")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        
        print("\nIf you're having issues, make sure you've run 'python generate_sample_data.py' first.")

if __name__ == "__main__":
    main() 
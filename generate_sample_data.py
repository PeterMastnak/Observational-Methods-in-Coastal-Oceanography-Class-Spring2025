#!/usr/bin/env python
"""
Generate a synthetic oceanographic data file in NetCDF format.
This will create a simple grid with temperature data for testing the process_ocean_data.py script.
"""

import numpy as np
import xarray as xr
import os
from datetime import datetime

# Create a data directory if it doesn't exist
data_dir = 'data'
os.makedirs(data_dir, exist_ok=True)

# Parameters for the grid
lon_min, lon_max = -80, -60  # Western Atlantic
lat_min, lat_max = 30, 50    # Mid-latitudes
depth_min, depth_max = 0, 1000
n_lon, n_lat, n_depth = 20, 20, 20

# Create coordinate arrays
lon = np.linspace(lon_min, lon_max, n_lon)
lat = np.linspace(lat_min, lat_max, n_lat)
depth = np.linspace(depth_min, depth_max, n_depth)
time = np.array([0])  # Single time point

# Create a synthetic temperature field
# Temperature decreases with depth and has some horizontal structure
temp = np.zeros((1, n_depth, n_lat, n_lon))

# Base temperature profile that decreases with depth
temp_profile = 25 * np.exp(-depth / 200) + 5

# Add horizontal structure
for i in range(n_lon):
    for j in range(n_lat):
        # Add some horizontal variation - warmer in the south, cooler in the north
        lat_factor = 1.0 - 0.5 * (lat[j] - lat_min) / (lat_max - lat_min)
        
        # Add a warm core eddy in the middle
        eddy_lon = -70  # Center longitude
        eddy_lat = 40   # Center latitude
        eddy_radius = 5 # Degrees
        dist = np.sqrt((lon[i] - eddy_lon)**2 + (lat[j] - eddy_lat)**2)
        eddy_factor = 5 * np.exp(-dist**2 / (2 * eddy_radius**2))
        
        # Combine factors
        for k in range(n_depth):
            temp[0, k, j, i] = temp_profile[k] * lat_factor + eddy_factor * np.exp(-depth[k] / 100)

# Create an xarray dataset
ds = xr.Dataset(
    data_vars={
        "t_an": (["time", "depth", "lat", "lon"], temp, {
            "long_name": "Temperature",
            "units": "degC",
            "standard_name": "sea_water_temperature",
        }),
    },
    coords={
        "lon": (["lon"], lon, {
            "long_name": "Longitude",
            "units": "degrees_east",
            "standard_name": "longitude",
        }),
        "lat": (["lat"], lat, {
            "long_name": "Latitude",
            "units": "degrees_north",
            "standard_name": "latitude",
        }),
        "depth": (["depth"], depth, {
            "long_name": "Depth",
            "units": "m",
            "standard_name": "depth",
            "positive": "down",
        }),
        "time": (["time"], time, {
            "long_name": "Time",
            "units": "days since 2000-01-01",
            "standard_name": "time",
        }),
    },
    attrs={
        "title": "Synthetic Oceanographic Data",
        "description": "Synthetic temperature data for testing oceanographic processing",
        "history": f"Created {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "source": "generate_sample_data.py",
    },
)

# Save to NetCDF file
output_file = os.path.join(data_dir, 'synthetic_ocean_temp.nc')
ds.to_netcdf(output_file)
print(f"Created synthetic ocean data file: {output_file}")
print(f"Dimensions: lon({n_lon}) x lat({n_lat}) x depth({n_depth})")
print(f"Temperature range: {temp.min():.2f} to {temp.max():.2f} °C") 
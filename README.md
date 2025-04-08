# Oceanographic Data Analysis Tools

This repository contains Python scripts for analyzing oceanographic data using the GSW-Python package (TEOS-10). It was developed as part of a graduate course in Observational Methods in Coastal Oceanography.

## Features

- Exploration of seawater density relationships with pressure, temperature, and salinity
- Analysis of water mass mixing and cabbeling effects
- Generation of synthetic oceanographic data
- Processing and visualization of ocean temperature and salinity profiles
- Comparison of different oceanographic measures (Absolute vs. Practical Salinity, Conservative vs. In-situ Temperature)

## Installation

1. Clone this repository:
```bash
git clone [repository-url]
```

2. Create and activate a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows, use: venv\Scripts\activate
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

## Scripts

- `explore_density_relationships.py`: Investigates how pressure, temperature, and salinity affect seawater density
- `gsw_examples.py`: Demonstrates basic usage of the GSW package and creates T-S diagrams
- `process_ocean_data.py`: Processes and visualizes oceanographic data
- `generate_sample_data.py`: Creates synthetic ocean data for testing and learning

## Usage

Each script can be run independently. For example:

```bash
python explore_density_relationships.py
```

This will generate various plots and a summary report in the current directory.

## Output Files

The scripts generate several visualization files:
- Density relationship plots (PNG format)
- Temperature and salinity profiles
- T-S diagrams
- Summary reports (TXT format)

## Requirements

- Python 3.7+
- gsw-python
- numpy
- matplotlib
- netCDF4
- cmocean

See `requirements.txt` for specific version requirements.

## Contributing

This is an educational project. Feel free to fork and modify for your own use. If you find any bugs or have suggestions, please open an issue.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Developed for the Observational Methods in Coastal Oceanography course
- Uses the GSW-Python package based on TEOS-10 
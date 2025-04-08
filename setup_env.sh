#!/bin/bash
# Script to set up the virtual environment for the Observational Methods in Coastal Oceanography project

echo "Setting up virtual environment for Coastal Oceanography project..."

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 could not be found. Please install Python 3 before continuing."
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
else
    echo "Virtual environment already exists."
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install required packages
echo "Installing required packages from requirements.txt..."
pip install -r requirements.txt

# Check if GSW installed correctly
echo "Verifying GSW installation..."
python -c "import gsw; print(f'GSW version: {gsw.__version__}')"

echo ""
echo "Setup complete! You can now run the example script using:"
echo "source venv/bin/activate"
echo "python gsw_examples.py"
echo ""
echo "To deactivate the virtual environment, simply run:"
echo "deactivate" 
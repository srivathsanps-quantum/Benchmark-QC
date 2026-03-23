#!/bin/bash
set -e  # Stop on any error

echo "Step 1: Installing pip packages..."
pip install -r requirements.txt

echo "Step 2: Cloning and installing ActiveSpaceFinder..."
git clone https://github.com/HQSquantumsimulations/ActiveSpaceFinder
cd ActiveSpaceFinder
pip install .
cd ..

echo "Step 3: Configuring dmrgscf settings..."
./ActiveSpaceFinder/init_dmrgscf_settings.sh

echo "All done! Environment is ready."
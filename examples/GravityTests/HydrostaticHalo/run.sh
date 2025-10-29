#!/bin/bash

# Allow passing arguments, otherwise use defaults
OUTPUT_DIR=${1:-results}
IC_FOLDER=${2:-None}

# Create output directory if it doesn't exist
if [ -d "$OUTPUT_DIR" ]; then
    echo "Output directory '$OUTPUT_DIR' already exists, using it."
else
    echo "Creating output directory '$OUTPUT_DIR'."
    mkdir -p "$OUTPUT_DIR"
fi

# --- Generate initial conditions if IC_FOLDER is defined ---
if [ "$IC_FOLDER" != "None" ] && [ -n "$IC_FOLDER" ]; then
    echo "Initial conditions folder defined: $IC_FOLDER"
    echo "Generating initial conditions using makeIC.py..."

    python3 makeIC.py 100000 "$OUTPUT_DIR" "$IC_FOLDER"

else
    # No IC folder specified
    if [ ! -e "$OUTPUT_DIR/Hydrostatic.hdf5" ]; then
        echo "No IC folder specified. Generating default initial conditions..."
        python3 makeIC.py 100000 "$OUTPUT_DIR" "$IC_FOLDER"
    else
        echo "Initial conditions already exist at $OUTPUT_DIR/Hydrostatic.hdf5."
    fi
fi

# Run SWIFT
echo "Running SWIFT, output will go to $OUTPUT_DIR"

# Get the directory of this script (absolute path, no symlinks)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Path to swift binary relative to script location
SWIFT_BIN="$SCRIPT_DIR/../../../swift"

# Path to YAML file (also relative to script location)
YAML_FILE="$SCRIPT_DIR/hydrostatic.yml"

# Run SWIFT from inside $OUTPUT_DIR, but with paths to ICs/YAML relative to parent
(
    cd "$OUTPUT_DIR"
    "$SWIFT_BIN" --external-gravity --hydro --threads=1 "$YAML_FILE" \
        2>&1 | tee output.log
)

echo "Plotting density profiles"
# make parent directory of plots
mkdir -p "$OUTPUT_DIR/plots/density_profile"
python3 density_profile.py 2. 200 300 "$OUTPUT_DIR"

echo "Plotting internal energy profiles"
mkdir "$OUTPUT_DIR/plots/internal_energy"
python3 internal_energy_profile.py 2. 200 300 "$OUTPUT_DIR"

echo "Plotting radial velocity profiles"
mkdir "$OUTPUT_DIR/plots/radial_velocity_profile"
python3 velocity_profile.py 2. 200 300 "$OUTPUT_DIR"

echo "Plotting energy as a function of time"
python3 test_energy_conservation.py 300 "$OUTPUT_DIR"

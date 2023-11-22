#!/bin/bash

# Batch job settings for UCL Legion
#$ -l h_rt=2:00:00
#$ -l h_vmem=2G
#$ -N uptrop_test_job_source
#$ -pe smp 1
#$ -V -cwd
# -P hpc.10
# -l paid=1

# Load environment modules
# export OMP_NUM_THREADS=1

# Python environment settings
CONDA_ENV_PATH="/Path/to/the/python/env/bin/activate"
PYTHONPATH_DIR="/Path/to/the/parent/directory/of/uptrop/"
export PYTHONPATH="$PYTHONPATH:$PYTHONPATH_DIR"

# Activate the virtual environment
source "$CONDA_ENV_PATH" || { echo "Failed to source conda environment"; exit 1; }

# Input and output directories
DATA_DIR="/lustre/projects/uptrop/tropomi/Data/"
OUTPUT_DIR="/Path/to/Output/Directory/"

# Date and pressure range settings
START_DATES=("2020-06-01")
END_DATES=("2020-08-31")
PMIN_VALUES=(450)
PMAX_VALUES=(600)

# Function to run the Python script
run_python_script() {
    species=$1
    product=$2
    start_date=$3
    end_date=$4
    pmin=$5
    pmax=$6

    python -m uptrop.main \
           --species "$species" \
           --product "$product" \
           --trop_dir "$DATA_DIR" \
           --out_dir "$OUTPUT_DIR" \
           --start_date "$start_date" \
           --end_date "$end_date" \
           --pmin $pmin \
           --pmax $pmax

    # Check if the Python script executed successfully
    if [ $? -ne 0 ]; then
        echo "Python script failed for species: $species, product: $product, start date: $start_date, end date: $end_date, PMIN: $pmin, PMAX: $pmax"
        exit 1
    fi
}

# Loop through date and pressure ranges and submit jobs
for i in "${!START_DATES[@]}"; do
    for j in "${!PMIN_VALUES[@]}"; do
        # Run the script for different species, products, date ranges, and pressure ranges
        run_python_script "no2" "PAL" "${START_DATES[$i]}" "${END_DATES[$i]}" "${PMIN_VALUES[$j]}" "${PMAX_VALUES[$j]}"
        run_python_script "o3" "OFFL" "${START_DATES[$i]}" "${END_DATES[$i]}" "${PMIN_VALUES[$j]}" "${PMAX_VALUES[$j]}"
    done
done

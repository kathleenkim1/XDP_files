#!/bin/bash
#SBATCH --job-name=map_cells         # Job name
#SBATCH --output=map_cells_%j.out    # Output log file (%j inserts job ID)
#SBATCH --error=map_cells_%j.err     # Error log file
#SBATCH --time=4:00:00               # Max runtime (HH:MM:SS)
#SBATCH --mem=128G                   # Memory (adjust as needed)
#SBATCH --cpus-per-task=16           # Number of CPU cores
#SBATCH --partition=disco          # Adjust to your cluster's partition

# Activate Conda environment (if needed)
source ~/.bashrc
conda activate cell_mapping_env  # Replace with your actual Conda env

# Run the mapping script
bash /broad/macosko/kimkathl/cell_type_mapper/run_mapmycells.sh --ref /broad/macosko/kimkathl/hmba_neuron.h5ad --query /broad/macosko/kimkathl/xdp_neurons.h5ad -o /broad/macosko/kimkathl/new_output

#!/bin/bash
#SBATCH --account <quest_allocation>
#SBATCH --partition short
#SBATCH --job-name maxatac_prepare
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 60GB
#SBATCH --time 2:00:00
#SBATCH --output logs/maxatac/prepare/%j.log
#SBATCH --verbose

# Load in modules
module purge
module load python-miniconda3/4.10.3

# Activate conda env
source activate maxatac

# Define dirs
output_dir="/path/to/output_dir"

# Run maxatac prepare
maxatac prepare -i "${1}" \
-o "${output_dir}" \
-prefix "${2}"
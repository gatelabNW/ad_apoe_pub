#!/bin/bash
#SBATCH --account <quest_allocation>
#SBATCH --partition normal
#SBATCH --job-name footprinting_batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 8:00:00
#SBATCH --output logs/footprinting/batch/%j.log
#SBATCH --verbose

# Describe arguments
# ${1} = cell type
# ${2} = motif

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load R/4.1.1 hdf5/1.8.15-serial geos/3.8.1

# Run job
Rscript code/footprinting/batch/footprinting-batch.R --celltype="${1}" --motif="${2}"

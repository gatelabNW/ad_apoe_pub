#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name da_batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 8:00:00
#SBATCH --output logs/da/batch/%j.log
#SBATCH --verbose

# Describe arguments
# ${1} = comparison
# ${2} = cell type

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load R/4.1.1 hdf5/1.8.15-serial geos/3.8.1

# Run job
Rscript code/da/batch/da-batch.R --comparison="${1}" --celltype="${2}"
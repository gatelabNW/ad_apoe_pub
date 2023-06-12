#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name DElegate
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 100GB
#SBATCH --time 1:00:00
#SBATCH --output logs/DElegate/batch/%j.log
#SBATCH --verbose

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load R/4.1.1 hdf5/1.8.15-serial geos/3.8.1

# Run job
Rscript code/DElegate/batch/DElegate-batch.R \
--comparison="${1}"
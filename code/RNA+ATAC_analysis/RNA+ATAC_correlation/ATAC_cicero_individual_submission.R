#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name cicero_batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 120GB
#SBATCH --time 12:00:00
#SBATCH --output logs/cicero/batch/%j.log
#SBATCH --verbose

# Describe arguments
# ${1} = diagnosis
# ${2} = geno
# ${3} = cell type

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load R/4.1.1 sqlite/3.27.2 proj/7.1.1 gdal/3.1.3-R-4.1.1 hdf5/1.8.15-serial geos/3.8.1

# Run job
Rscript code/cicero/batch/cicero-batch.R --diagnosis="${1}" --geno="${2}" --celltype="${3}"
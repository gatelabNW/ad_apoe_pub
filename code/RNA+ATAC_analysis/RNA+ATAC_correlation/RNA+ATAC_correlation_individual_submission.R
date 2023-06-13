#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name exp_acc_corr_batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 2:00:00
#SBATCH --output logs/expression_accessibility_correlation/batch/%j.log
#SBATCH --verbose

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load R/4.1.1 sqlite/3.27.2 proj/7.1.1 gdal/3.1.3-R-4.1.1 hdf5/1.8.15-serial geos/3.8.1

# Run job
Rscript code/expression_accessibility_correlation/batch/expression_accessibility_correlation-batch.R \
--atac_celltype="${1}" \
--rna_celltype="${2}" \
--coacc_cutoff="${3}" \
--diagnosis="${4}"
#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name maxatac_peaks
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 60GB
#SBATCH --time 2:00:00
#SBATCH --output logs/maxatac/peaks/%j.log
#SBATCH --verbose

# ${1} = tf
# ${2} = bigwig
# ${3} = prefix
# ${4} = output_dir

# Print date and name
date
echo 'Natalie Piehl'

# Load in modules
module purge
module load python-miniconda3/4.10.3

# Activate conda env
source activate maxatac

# Run maxatac prepare
maxatac peaks \
-i "${2}" \
-cutoff_file "~/opt/maxatac/data/models/${1}/${1}_validationPerformance_vs_thresholdCalibration.tsv" \
--output "${4}" \
--prefix "${3}_${1}"

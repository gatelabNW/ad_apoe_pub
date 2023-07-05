#!/bin/bash
#SBATCH --account <quest_allocation>
#SBATCH --partition short
#SBATCH --job-name maxatac_tf_enrichment
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 180GB
#SBATCH --time 4:00:00
#SBATCH --output logs/maxatac/tf_enrichment/%j.log
#SBATCH --verbose

# ${1} = tf
# ${2} = bigwig
# ${3} = prefix
# Load in modules
module purge
module load python-miniconda3/4.10.3

# Activate conda env
source activate maxatac

# Define output dir
output_dir="/path/to/output_dir/"
[ -d "$output_dir" ] || mkdir -p "$output_dir"

echo "Processing ${1} in ${3}"

# Run maxatac predict
maxatac predict \
-tf "${1}" \
--signal "${2}" \
-o "${output_dir}" \
--prefix "${3}_${1}"
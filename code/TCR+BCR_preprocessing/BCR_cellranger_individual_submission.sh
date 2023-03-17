#!/bin/bash
#SBATCH --account <gate_lab_allocation>
#SBATCH --partition normal
#SBATCH --job-name cellranger_bcr
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 08:00:00
#SBATCH --output /path/to/logs/cellranger/bcr/%x-%j.log
#SBATCH --verbose

# ${1} = id
# ${2} = sample_dir
# ${3} = sample
# ${4} = output_dir

date

# Module prep
module purge
module load cellranger/6.0.0

# Navigate to output directory
cd ${4}

# Run cellranger vdj
cellranger vdj \
--id ${1} \
--fastqs ${2} \
--chain "IG" \
--reference "/path/to/ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \
--sample ${3} \
--localcores 16

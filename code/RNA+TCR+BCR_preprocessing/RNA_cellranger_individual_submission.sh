#!/bin/bash
#SBATCH --account <gate_lab_allocation>
#SBATCH --partition genomics
#SBATCH --job-name cellranger_count
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 8:00:00
#SBATCH --output /path/to/logs/cellranger/rna/%x-%j.log
#SBATCH --verbose

# ${1} = id
# ${2} = sample_dir
# ${3} = sample
# ${4} = cell_num
# ${5} = output_dir

date

# Module prep
module purge
module load cellranger/6.0.0

# Navigate to output directory
cd ${5}

# Run cellranger count
cellranger count \
--id ${1} \
--fastqs ${2} \
--transcriptome "/path/to/ref/refdata-gex-GRCh38-2020-A" \
--sample ${3} \
--expect-cells ${4} \
--localcores 8

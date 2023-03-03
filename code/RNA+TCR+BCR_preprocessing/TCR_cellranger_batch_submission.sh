# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-28-2022
# Written by: Natalie Piehl
# Summary: Run cellranger on scTCR fastqs
#
#-------------------------------------------------------------------------------

# Define and check directory containing GEX fastq files
fastq_dir="/projects/b1042/Gate_Lab/fromNUSeq/Gate02_9.22.2022/"
[ -d "$fastq_dir" ] && echo "Directory to fastq files exists."

# Define, make (if necessary), and navigate to output directory
output_dir="/projects/b1042/Gate_Lab/AD_APOE/results/cellranger/tcr/"
[ -d "$output_dir" ] || mkdir "$output_dir"

# Create array of fastq_dir subdirectories
sample_arr=()
while IFS=  read -r -d $'\0'; do
   sample_arr+=("$REPLY")
done < <(find "$fastq_dir" -maxdepth 1 -mindepth 1 -type d -print0)

# For each sample in fastq directory...
# for sample_dir in "${sample_arr[0]}"; do
for sample_dir in "${sample_arr[@]:1}"; do
    # Isolate sample ID
    id=$(basename "$sample_dir")
    echo "Submitting $id now"

    # Run cellranger vdj script
    sbatch batch_cellranger_tcr.sh "$id" "$sample_dir" "$id" "$output_dir"
    sleep 2
done

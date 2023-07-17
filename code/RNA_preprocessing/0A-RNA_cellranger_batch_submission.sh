# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 07-25-2022
# Written by: Natalie Piehl
# Summary: Run cellranger on scRNA fastqs
#
#-------------------------------------------------------------------------------

# Define and check directory containing GEX fastq files
fastq_dir="/path/to/fastqs"
[ -d "$fastq_dir" ] && echo "Directory to fastq files exists."

# Define, make (if necessary), and navigate to output directory
output_dir="/path/to/output/folder"
[ -d "$output_dir" ] || mkdir "$output_dir"

# Create array of fastq_dir subdirectories
sample_arr=()
while IFS=  read -r -d $'\0'; do
    sample_arr+=("$REPLY")
done < <(find "$fastq_dir" -maxdepth 1 -mindepth 1 -type d -print0)

counter=0
# For each sample in fastq directory...
for sample_dir in "${sample_arr[@]}"
# for sample_dir in "${sample_arr[0]}"
do
    let counter++
    # if [ $counter -ne 31 ]; then
    #   continue
    # fi

    # Isolate sample ID
    id=$(basename "$sample_dir")

    # Check if output directory already exists, and if so skip to next sample
    if [ -d "$output_dir/$id" ]; then
      continue
    fi

    # Get names of files in sample_dir
    files=()
    while IFS=  read -r -d $'\0'; do
        files+=("$REPLY")
    done < <(find "$sample_dir" -type f -name "*.fastq.gz" -print0)

    # Identify if need to cutoff last 24 or 25 characters
    sample_num=$(echo "$files[0]" | grep -o -P '(?<=_S).*(?=_L)')
    if [ ${#sample_num} -eq 2 ]; then
      cutoff=25; else
      cutoff=24
    fi

    # Isolate sample name
    sample=${files[0]::-$cutoff}
    sample=$(basename "$sample")
    echo "Submitting #$counter $id: $sample now"

    # Get predicted cell number from processing document
    cell_num=10000

    # Run cellranger count script
    sbatch batch_cellranger-rna.sh "$id" "$sample_dir" "$sample" "$cell_num" "$output_dir"
    sleep 2
done

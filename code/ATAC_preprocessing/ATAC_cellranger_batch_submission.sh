# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-26-2022
# Written by: Natalie Piehl
# Summary: Run cellranger atac on scATAC fastqs
#
#-------------------------------------------------------------------------------

# Specify paths
fastq_dir="/path/to/fastqs"
output_dir="/path/to/output/folder"
[ -d "$output_dir" ] || mkdir -p "$output_dir"

# Create array of fastq_dir subdirectories
sample_arr=()
while IFS=  read -r -d $'\0'; do
    sample_arr+=("$REPLY")
done < <(find "$fastq_dir" -maxdepth 1 -mindepth 1 -type d -print0)

counter=0
# Submit job for each sample
for sample_dir in "${sample_arr[@]}"
# for sample_dir in "${sample_arr[0]}"
do
  # Isolate sample ID
  id=$(basename "$sample_dir")

  let counter++
  # if [ $counter -ne 3 ]; then
  #   continue
  # fi

  # Check if output directory already exists, and if so skip to next sample
  if [ $counter -ne 17 ]; then
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

  # Run alignment job
  sbatch "batch_atac.sh" "$sample_dir" "$sample" "$id" "$output_dir"
  sleep 2
done

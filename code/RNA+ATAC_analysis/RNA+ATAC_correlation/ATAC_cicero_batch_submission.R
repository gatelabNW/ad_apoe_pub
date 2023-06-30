#!/bin/bash

# Define path to cell type color map
cell_type_colors_path="/path/to/celltype_color_map"

# Create array of cell types from metadata
cell_types=()
while IFS=',' read -ra row; do
  cell_types+=("${row[0]// /_}")
done < $cell_type_colors_path
cell_types=("${cell_types[@]:1}")

# Specify comparisons
diagnoses=("all" "Healthy Control" "Alzheimers Disease")
genos=("all" "E3/E3" "E3/E4" "E4/E4")

for diagnosis in "${diagnoses[@]}"
do
  for geno in "${genos[@]}"
  do
    # Process low mem cell types
    for cell_type in "${cell_types[0]}"
    do        
        # Print comparison and cell type
        cell_type="CD8+_T_Cells"
        echo "Submitting $diagnosis $geno $cell_type"

        # Submit individual job
        sbatch code/batch_jobs/cicero_individual_submission.sh "$diagnosis" "$geno" "$cell_type"

        # Sleep to give scheduler a break
        sleep 3
    done
  done
done

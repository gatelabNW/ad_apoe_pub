#!/bin/bash

# Define path to cell type color map
cell_type_colors_path="/path/to/celltype_color_map.csv"

# Create array of cell types from metadata
cell_types=()
while IFS=',' read -ra row; do
  cell_types+=("${row[0]// /_}")
done < $cell_type_colors_path
cell_types=("${cell_types[@]:1}")

# Specify long time requiring cell types
long_cell_types=("CD4+_T_Cells" "NK_Cells" "CD8+_T_Cells")

# Remove long time cell types from list
for i in "${long_cell_types[@]}"; do
  cell_types=(${cell_types[@]//*$i*})
done

# Specify comparisons
comparisons=("ADvsHC" \
"ADvsHC_33" "ADvsHC_34" "ADvsHC_44" \
"34vs33_HC" "34vs33_AD" \
"44vs33_HC" "44vs33_AD" \
"44vs34_HC" "44vs34_AD" \
)

for comparison in "${comparisons[9]}"
do
    # Process low mem cell types
    for cell_type in "${cell_types[@]}"
    do        
        # Print comparison and cell type
        echo "Submitting $comparison in $cell_type for 8 hours"

        # Submit individual job
        sbatch code/batch_jobs/da_individual_submission.sh "$comparison" "$cell_type"

        # Sleep to give scheduler a break
        sleep 3
    done

    # Process long time cell types
    for cell_type in "${long_cell_types[@]}"
    do        
        # Print comparison and cell type
        echo "Submitting $comparison in $cell_type for 24 hours"

        # Submit individual job
        sbatch code/batch_jobs/da_individual_submission_long.sh "$comparison" "$cell_type"

        # Sleep to give scheduler a break
        sleep 3
    done
done
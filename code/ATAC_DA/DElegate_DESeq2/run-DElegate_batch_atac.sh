#!/bin/bash

# Define comparisons
comparisons=("ADvsHC" \
"ADvsHC_33" "ADvsHC_34" "ADvsHC_44" \
"ADvsHC_4carrier" "4carriervs33_AD" "4carriervs33_HC" \
"34vs33_HC" "34vs33_AD" \
"44vs33_HC" "44vs33_AD" \
"44vs34_HC" "44vs34_AD" \
"34vs33" "44vs33" "44vs34" "4carriervs33"
)

for comparison in "${comparisons[@]}"
do
      # Print comparison and cell type
      echo "Submitting $comparison"

      # Run job
      sbatch code/DElegate/batch_atac/DElegate_batch_atac.sh "$comparison"

      # Sleep to give scheduler a break
      sleep 3
done
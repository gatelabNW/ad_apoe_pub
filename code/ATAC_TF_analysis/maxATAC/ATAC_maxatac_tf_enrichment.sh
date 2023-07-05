#!/bin/bash

# Print date and name
date
echo 'Natalie Piehl'

# Define dirs
input_dir="/path/to/maxatac/average/results/"

# Specify cell types
celltypes=("B_Cells" "Monocytes" "Dendritic_Cells" \
            "CD4+_T_Cells" "CD8+_T_Cells" "NK_Cells")

# Identify all TFs
tf_dir="/home/ncp2306/opt/maxatac/data/models/"
tf_paths=()
while IFS=  read -r -d $'\0'; do
    tf_paths+=("$REPLY")
done < <(find "$tf_dir" -maxdepth 1 -mindepth 1 -type d -print0)

for celltype in "${celltypes[4]}"
do
    # Identify cell type specific bigwigs
    bigwigs=()
    while IFS=  read -r -d $'\0'; do
        bigwigs+=("$REPLY")
    done < <(find "$input_dir" -maxdepth 1 -mindepth 1 -type f -name "*${celltype}*" -print0)

    # Iterate through bigwig files
    for bigwig in "${bigwigs[@]:3:1}"
    do
        # Isolate prefix from bigwig
        base=$(basename $bigwig)
        prefix=${base%".bw"}
        echo $prefix

        for tf_path in "${tf_paths[@]}"
        do
            # Print submission
            tf=$(basename $tf_path)
            echo "Submitting $tf in $prefix"

            # Submit job
            sbatch /projects/b1169/nat/AD_APOE/code/maxatac/tf_enrichment/maxatac-all_tf_enrichment.sh \
            "$tf" "$bigwig" "$prefix"

            # Sleep to give scheduler a break
            sleep 3
        done
    done
done



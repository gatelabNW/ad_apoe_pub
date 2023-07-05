#!/bin/bash

# Define dirs
input_dir="/path/to/tf_enrichment/results/"
output_dir="/path/to/output_dir/"
[ -d "$output_dir" ] || mkdir -p "$output_dir"

# Specify cell types
celltypes=("B_Cells" "Monocytes" "Dendritic_Cells" \
            "CD4+_T_Cells" "CD8+_T_Cells" "NK_Cells")

for celltype in "${celltypes[@]}"
do
    # Identify cell type specific bigwigs
    bigwigs=()
    while IFS=  read -r -d $'\0'; do
        bigwigs+=("$REPLY")
    done < <(find "$input_dir" -maxdepth 1 -mindepth 1 -type f -name "*${celltype}*.bw" -print0)

    # Iterate through bigwig files
    for bigwig in "${bigwigs[@]:0:700}"
    do
        # Isolate prefix and tf from bigwig
        base=$(basename $bigwig)
        if [[ "$bigwig" == *"4carrier"* ]];then
            prefix=$(echo $base | cut -d_ -f -3)
            tf=$(echo $base | cut -d_ -f4)
        else
            prefix=$(echo $base | cut -d_ -f -2)
            tf=$(echo $base | cut -d_ -f3)
        fi
        tf=${tf%.bw}
        

        # Check if output already exists
        if [ -f "${output_dir}${prefix}_${tf}_32bp.bed" ];then
            continue
        fi

        echo "Submitting $tf in $prefix"

        # Submit job
        sbatch /projects/b1169/nat/AD_APOE/code/maxatac/peaks/maxatac-peaks.sh \
        "$tf" "$bigwig" "$prefix" "$output_dir"

        # Sleep to give scheduler a break
        sleep 2
    done
done



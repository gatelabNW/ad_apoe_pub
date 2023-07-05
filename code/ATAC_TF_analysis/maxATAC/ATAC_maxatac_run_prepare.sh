#!/bin/bash

# Define input dir
input_dir="/path/to/pseudobulk/fragments"

# Generate array of pseudobulk fragments
frag_arr=()
while IFS=  read -r -d $'\0'; do
    frag_arr+=("$REPLY")
done < <(find "$input_dir" -maxdepth 1 -mindepth 1 -type f -print0)

# 1-20 submitted
counter=0
for frag in "${frag_arr[@]}"
do
    # Update counter
    let counter++
    if (($counter >= 1 && $counter <= 100)); then
        # Define prefix
        base=$(basename $frag)
        prefix=${base%"_fragments.tsv"}

        echo "$counter: $prefix"

        # Run maxatac prepare
        sbatch code/maxatac/prepare/maxatac-ind_prepare.sh "${frag}" "${prefix}"

        # Sleep for scheduler
        sleep 3
    fi
done
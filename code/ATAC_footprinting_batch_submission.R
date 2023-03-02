#!/bin/bash

# Specify cell type
cell_type="NK_Cells"

# Specify motifs
motifs=("TCFL5" "NRF1" "ZBTB14" "ZBTB33" "CTCFL" "Wt1" "Klf1" "CTCF" "FOSB::JUNB" "RUNX2" "FOS" "MTF1" "Zic2" \
"FOS::JUNB" "FOS::JUND" "RUNX1" "RELA" "REL" "ETV1" "GABPA" "ELF1" "ETV4" "RREB1" "RUNX3")

for motif in "${motifs[@]}"
do
    # Print motif
    echo "Submitting $motif in $cell_type"

    # Submit individual job
    sbatch code/batch_jobs/footprinting_individual_submission.sh "$cell_type" "$motif"

    # Sleep to give scheduler a break
    sleep 3
done
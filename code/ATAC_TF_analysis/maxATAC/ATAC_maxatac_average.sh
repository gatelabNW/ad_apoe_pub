#!/bin/bash
#SBATCH --account <quest_allocation>
#SBATCH --partition normal
#SBATCH --job-name maxatac_average
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 100GB
#SBATCH --time 8:00:00
#SBATCH --output logs/maxatac/average/%j.log
#SBATCH --verbose

# Load in modules
module purge
module load python-miniconda3/4.10.3

# Activate conda env
source activate maxatac

# Define dirs
input_dir="/path/to/maxatac/prepared/data/"
output_dir="/path/to/output_dir/"

# Specify diagnosis + apoe combination
diagnosis_apoe="AD"

# Select for appropriate samples
if [ $diagnosis_apoe == "AD" ]; then
    samples=("A773_y4" "A1055_y2" "A1092_y2" "A230_y3" "A1154" "A942" "A516_y5" "A906_Y3" \
    "A1052_y4" "A921_y3" "A1120_y2" "A932_y3" "A1160" "A917_y2" "A656_y3" "A1147" "A863_y2" "A1241_Y3" "A1236" \
    "A1034" "A254_y4" "A70" "A802_y4" "A947" "A965_y3" "A738")
elif [ $diagnosis_apoe == "HC" ]; then
    samples=("A820_y4" "A659_y2" "A598_y2" "A968_y2" "A912" "A1086_Y2" "A836" "A1028" \
    "A1020_y2" "A905_y2" "A1162" "A781_y4" "A1282" "A911_y2" "A1200" "A1180" "A1081_Y2" \
    "A960_y2" "A978_y2" "A780_y4" "A989" "A970_y2" "A1279" "A1010_y2")
fi

# Generate array of prepared pseudobulks
frag_arr=()
while IFS=  read -r -d $'\0'; do
    frag_arr+=("$REPLY")
done < <(find "$input_dir" -maxdepth 1 -mindepth 1 -type f -name "*_minmax01.bw" -print0)

# Specify cell types
celltypes=("B_Cells" "Monocytes" "Dendritic_Cells" \
            "CD4+_T_Cells" "CD8+_T_Cells" "NK_Cells" \
            "Other_T_Cells" "Other")

# Iterate through cell types
for celltype in "${celltypes[@]}"
do
    # # Identify files to merge
    sample_files=()
    for file in "${frag_arr[@]}"; do
        for sample in "${samples[@]}"; do
            if [[ "$file" == *"$sample"* ]] && [[ "$file" == *"$celltype"* ]]; then
                sample_files+=("$file")
            fi 
        done
    done

    # Run maxatac prepare
    maxatac average \
    -i "${sample_files[@]}" \
    --output "${output_dir}" \
    --prefix "${diagnosis_apoe}_${celltype}"
done

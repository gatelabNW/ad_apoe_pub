#!/bin/bash

# List parameters
# 1:0-2, 2:3-5, 3:6-7
atac_cell_types=("B_Cells" "Monocytes" "Dendritic_Cells" \
                 "CD4+_T_Cells" "CD8+_T_Cells" "NK_Cells" \
                 "Other_T_Cells" "Other")
# 1:0-3, 2:4-9, 3:10-15, 4:16-20, 5:21-23, 6:24-26, 7:27-30
rna_cell_types=("B_intermediate" "B_memory" "B_naive" "Plasmablast" \
                "CD14_Mono" "CD16_Mono" "ASDC" "cDC1" "cDC2" "pDC" \
                "CD4_CTL" "CD4_Naive" "CD4_Proliferating" "CD4_TCM" "CD4_TEM" "Treg" \
                "CD8_Naive" "CD8_Proliferating" "CD8_TCM" "CD8_TEM" "MAIT" \
                "NK" "NK_Proliferating" "NK_CD56bright" \
                "ILC" "dnT" "gdT" \
                "Platelet" "Eryth" "HSPC" "Doublet")
coacc_cutoff=0.01                
diagnoses=("Healthy_Control" "Alzheimers_Disease")
genos=("E3/E3" "E3/E4" "E4/E4")

for atac_cell_type in "${atac_cell_types[@]}"
do
    if [[ $atac_cell_type == "B_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:0:4}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "Monocytes" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:4:2}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "Dendritic_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:6:4}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "CD4+_T_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:10:6}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "CD8+_T_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:16:5}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "NK_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:21:3}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "Other_T_Cells" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:24:3}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi

    if [[ $atac_cell_type == "Other" ]];
    then
        for rna_cell_type in "${rna_cell_types[@]:27:3}"
        do
            for diagnosis in "${diagnoses[@]}"
            do
                # Print comparison and cell type
                echo "Submitting $diagnosis $atac_cell_type $rna_cell_type"

                # Run job
                sbatch code/batch_jobs/exp_acc_corr_individual_diagnosis_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff" "$diagnosis"

                # Sleep to give scheduler a break
                sleep 3
            done
            
            # Print comparison and cell type
            echo "Submitting HC and AD in $atac_cell_type $rna_cell_type"

            # Run job
            sbatch code/batch_jobs/exp_acc_corr_individual_submission.sh "$atac_cell_type" "$rna_cell_type" "$coacc_cutoff"

            # Sleep to give scheduler a break
            sleep 3

        done
    fi
done
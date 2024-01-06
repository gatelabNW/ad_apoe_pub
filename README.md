# Epigenetic dysregulation in Alzheimer’s disease peripheral immunity
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

![updated_schematic](https://github.com/gatelabNW/ad_apoe_pub/assets/91904251/8f69bc36-3fc5-44da-b49b-b3c8e395b438)

## Abstract

Much attention has been paid to the role of the perenchymal brain immune response in Alzheimer’s disease (AD). Yet, the peripheral immune system in AD has not been thoroughly studied with modern sequencing methods. Here, we used a combination of single-cell sequencing strategies, including assay for transposase-accessible chromatin and RNA sequencing, to investigate the epigenetic and transcriptional alterations to the AD peripheral immune system. We reveal a striking amount of open chromatin in peripheral immune cells in AD. In CD8 T cells, we uncover a cis-regulatory DNA element co-accessible with the CXC motif chemokine receptor 3 gene promoter. In monocytes, we identify a novel RELA transcription factor binding site adjacent to an open chromatin region in the nuclear factor kappa B subunit 2 gene. We also demonstrate apolipoprotein E genotype-dependent epigenetic changes that correspond to altered cytokine gene expression in monocytes. Surprisingly, we also identify differentially accessible chromatin regions in genes associated with sporadic AD. Finally, we provide an online data portal resource to explore gene expression in the AD immune system. Our findings provide novel insights into the complex relationship between epigenetics and genetic risk factors in the AD peripheral immune system.

## About
This repository contains code used to process and analyze scRNA+TCR+BCRseq and scATACseq data from the **Epigenetic dysregulation in Alzheimer’s disease peripheral immunity** study. Preprocessing scripts are numbered to specify execution order and corresponding scripts for figure panels are listed below.

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```ad_apoe_pub.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The RNA+BCR+TCR dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/ad_apoe_rna/).

Raw RNA+BCR+TCR .fastq files and gene expression matrices can be downloaded from [GSE226602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226602). Raw ATAC .fastq files and peak matrices can be downloaded from [GSE226267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226267).

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

## Figure Outline
All scripts in `code` folder.

#### Fig1: The peripheral immune system has more open chromatin in AD.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_preprocessing/2-ATAC_archr_umap.R` + `RNA_preprocessing/2-RNA_seurat_preprocessing.R`

#### Fig2: Epigenetic dysregulation in AD peripheral immunity and concordance with differential gene expression.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `ATAC_DA/ATAC_LR+DESeq2_overlap_barplot.R` + `ATAC_DA/ATAC_LR+DESeq2_overlap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_preprocessing/ATAC_DAR_peak_type_composition.R` + `ATAC_preprocessing/ATAC_DAR_up_down_composition.R`         
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA_DE/RNA_MAST+edgeR_overlap_barplot.R` + `RNA_DE/RNA_MAST+edgeR_overlap.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_barplot_upset.R`   

#### Fig3: Epigenetic changes to NF-κB signaling molecules in peripheral AD monocytes. 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_TF_analysis/ATAC_TF_enrichment_scatter.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C-D: `ATAC_TF_analysis/maxATAC/ATAC_ADvsHC_Monocytes_NFKB2.R`    

#### Fig4: Cis-regulatory sequences correlate with quantity of transcription of select CD8 T cell and monocyte genes.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_correlation/ATAC_cicero_summary.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `RNA+ATAC_analysis/RNA+ATAC_correlation/RNA+ATAC_cre-linked-gene+DEG_overlap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `RNA+ATAC_analysis/RNA+ATAC_correlation/RNA+ATAC_CXCR3_correlation.R`  

#### Fig5: APOE genotype-dependent innate immune dysregulation in AD.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_ADvsHC_in_APOE_heatmap_upset.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_ADvsHC_heatmap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `Misc/ATAC_coverage_plot.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `ATAC_TF_analysis/ATAC_TF_enrichment_scatter.R`   

#### Fig6: Epigenetic dysregulation of AD risk genes in the peripheral immune system.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A-D: `ATAC_DA/ATAC_DAR_AD_risk_gene_overlap.R` (coverage plots merged with heatmaps in Illustrator)      
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `Misc/ATAC_coverage_plot.R`     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F: `TCR_analysis/umap-clonal.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G: `RNA_DE/volcano.R`   

#### SuppFig1: Quality control metrics for scATACseq data.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `Misc/demographics/demographic_plots_optimized.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `Misc/revisions/biomarker_boxplots.R`
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C-E: `ATAC_preprocessing/1-ATAC_archr_preprocessing.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F: `ATAC_preprocessing/3-ATAC_archr_label_transfer.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G: `ATAC_preprocessing/1-ATAC_archr_preprocessing.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;H-J: `Misc/revisions/neuron_figures.R`

#### SuppFig2: Quality control metrics for scATACseq data.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA_preprocessing/4-RNA_celltype_annotation.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B-D: `RNA_preprocessing/2-RNA_seurat_preprocessing.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `RNA_preprocessing/1-RNA_soupx.R` (visual made with Prism)   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F: `RNA_preprocessing/5-seurat_rna_celltypemarkers_heatmap.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G: `RNA_preprocessing/3-RNA_integration_dim_reduc.R` + `Misc/revisions/neuron_figures.R` 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;I: `Misc/revisions/figures-cell_num_boxplots.R`

#### SuppFig3: Transcription factor footprinting and TFBS analysis in monocytes.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `ATAC_TF_analysis/ATAC_footprinting.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_TF_analysis/maxATAC/ATAC_ADvsHC_Monocytes_32bp_barplot_wholegenome.R`  

#### SuppFig5: 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_correlation/RNA+ATAC_cre-linked-gene+DEG_overlap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `RNA+ATAC_analysis/RNA+ATAC_correlation/RNA+ATAC_ABCA1_correlation.R`  

#### SuppFig7: APOE genotype-dependent dysregulation of clonally expanded T cells in AD.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `ATAC_DA/ATAC_DAR_AD_risk_gene_overlap.R` (coverage plots merged with heatmaps in Illustrator)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `TCR_analysis/tcr-clonal_DE_ADvsHC_within_APOE_genotype.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `Misc/ATAC_coverage_plot.R`  

#### SuppFig8: Clonally expanded CD8+ TCM cells are transcriptionally altered in AD.
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA_DE/RNA_DEG_APOE_heatmap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `RNA_DE/volcano.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `TCR_analysis/tcr-clonal_DE_APOEvsAPOE.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA_DE/volcano.R`  


![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/221924737-8ff64f66-bc81-4155-94a3-05121b393bfc.png)


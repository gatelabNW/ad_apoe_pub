# Apolipoprotein E genotype-dependent epigenetic immune dysregulation in Alzheimer’s disease
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

![merged_about_fig](https://user-images.githubusercontent.com/91904251/221924963-917c5637-83b7-4087-88f3-8316afdba8a7.png)

## Abstract

The apolipoprotein E (*APOE*) gene is a well-established risk factor for Alzheimer’s disease (AD). Much attention has been paid to the role of *APOE* in the brain immune response in AD. Yet, the relationship between various *APOE* genotypes and peripheral immunity in AD has not been studied.  Here, we used a combination of single cell sequencing strategies, including assay for transposase-accessible chromatin and RNA sequencing, to investigate the epigenetic and transcriptional influence of *APOE* genotypes on the AD immune system. We reveal *APOE* allele-dependent epigenetic changes that correspond to altered gene expression in peripheral monocytes and memory CD8 T cells in AD. We also identify differentially accessible chromatin regions in AD risk genes in peripheral immune cells. Finally, we provide an online data portal to explore gene expression in the AD immune system. Our findings provide novel insights into the complex relationship between *APOE* and the immune system in AD, and suggest potential therapeutic targets for this devastating disorder.

## About
This repository contains code used to process and analyze scRNA+TCR+BCRseq and scATACseq data from the **Peripheral immunomodulatory effects of Apolipoprotein E4 in Alzheimer’s disease** study. Preprocessing cripts are numbered to specify execution order and corresponding scripts for figure panels are listed below.

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```ad_apoe_pub.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The RNA+BCR+TCR dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/ad_apoe_rna/).

Raw RNA+BCR+TCR .fastq files and gene expression matrices can be downloaded from [GSE226602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226602). Raw ATAC .fastq files and peak matrices can be downloaded from [GSE226267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226267).

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

## Figure Outline
All scripts in `code` folder.

### Fig1: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_preprocessing/2-ATAC_archr_umap.R`  

### Fig2: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `ATAC_DA/ATAC_LR+DESeq2_overlap_barplot.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_preprocessing/ATAC_DAR_peak_type_composition.R` + `ATAC_preprocessing/ATAC_DAR_up_down_composition.R`         
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA_DE/RNA_MAST+edgeR_overlap_barplot.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `RNA_DE/RNA_MAST+edgeR_overlap.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_barplot_upset.R`   

### Fig3: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_TF_analysis/ATAC_TF_enrichment_scatter.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C+D: `ATAC_TF_analysis/maxATAC/ATAC_ADvsHC_Monocytes_NFKB2.R`    

### Fig4: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: Cicero promoter connections heatmap    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: Sig connections heatmap  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: cre-linked gene + DEG overlap cd14 monocytes  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: cre-linked gene + DEG overlap cd8 tems  

### Fig5: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A+B: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_ADvsHC_in_APOE_heatmap_upset.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_ADvsHC_heatmap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E+F: Chromatin tracks  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G: `ATAC_TF_analysis/ATAC_TF_enrichment_scatter.R`   

### Fig6: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA+ATAC_analysis/RNA+ATAC_DAR+DEG_overlap_scatter.R`       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: Chromatin tracks  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: AD risk gene heatmap   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: Tracks and heatmaps for AD risk genes   

### SuppFig1: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `demographics/demographic_plots_optimized.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B-E: `ATAC_preprocessing/1-ATAC_archr_preprocessing.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F: `ATAC_preprocessing/3-ATAC_archr_label_transfer.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;G: `ATAC_preprocessing/1-ATAC_archr_preprocessing.R`    

### SuppFig2: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA_preprocessing/4-RNA_celltype_annotation.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B-D: `RNA_preprocessing/2-RNA_seurat_preprocessing.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E: `RNA_preprocessing/1-RNA_soupx.R` (visual made with Prism)   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;F+G: `RNA_preprocessing/4-RNA_celltype_annotation.R`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;H: `ATAC_preprocessing/2-ATAC_archr_umap.R`    

### SuppFig3: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `ATAC_TF_analysis/ATAC_footprinting.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `ATAC_TF_analysis/maxATAC/ATAC_ADvsHC_Monocytes_32bp_barplot_wholegenome.R`  

### SuppFig4: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: ADvsHC DAR AD risk gene heatmap  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `TCR_analysis/umap-clonal.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `TCR_analysis/tcr-clonal_DE_ADvsHC_within_APOE_genotype.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA_DE/volcano.R`  

### SuppFig5: _______
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A: `RNA_DE/RNA_DEG_APOE_heatmap.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;B: `RNA_DE/volcano.R`  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C: `TCR_analysis/tcr-clonal_DE_APOEvsAPOE.R`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;D: `RNA_DE/volcano.R`  

![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/221924737-8ff64f66-bc81-4155-94a3-05121b393bfc.png)


# Apolipoprotein E genotype-dependent epigenetic immune dysregulation in Alzheimer’s disease
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

![merged_about_fig](https://user-images.githubusercontent.com/91904251/221924963-917c5637-83b7-4087-88f3-8316afdba8a7.png)

## Abstract

The apolipoprotein E (*APOE*) gene is a well-established risk factor for Alzheimer’s disease (AD). Much attention has been paid to the role of *APOE* in the brain immune response in AD. Yet, the relationship between various *APOE* genotypes and peripheral immunity in AD has not been studied.  Here, we used a combination of single cell sequencing strategies, including assay for transposase-accessible chromatin and RNA sequencing, to investigate the epigenetic and transcriptional influence of *APOE* genotypes on the AD immune system. We reveal *APOE* allele-dependent epigenetic changes that correspond to altered gene expression in peripheral monocytes and memory CD8 T cells in AD. We also identify differentially accessible chromatin regions in AD risk genes in peripheral immune cells. Finally, we provide an online data portal to explore gene expression in the AD immune system. Our findings provide novel insights into the complex relationship between *APOE* and the immune system in AD, and suggest potential therapeutic targets for this devastating disorder.

## About
This repository contains code used to process and analyze scRNA+TCR+BCRseq and scATACseq data from the **Peripheral immunomodulatory effects of Apolipoprotein E4 in Alzheimer’s disease** study. 

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```ad_apoe_pub.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The RNA+BCR+TCR dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/ad_apoe_rna/).

Raw RNA+BCR+TCR .fastq files and gene expression matrices can be downloaded from [GSE226602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226602). Raw ATAC .fastq files and peak matrices can be downloaded from [GSE226267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226267).

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/221924737-8ff64f66-bc81-4155-94a3-05121b393bfc.png)


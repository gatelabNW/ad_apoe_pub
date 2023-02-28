# Peripheral immunomodulatory effects of Apolipoprotein E4 in Alzheimer’s disease
[Gate Lab at Northwestern University](https://sites.northwestern.edu/gatelab/)

![merged_about_fig](https://user-images.githubusercontent.com/91904251/221924963-917c5637-83b7-4087-88f3-8316afdba8a7.png)

## Abstract

**Background:** Apolipoprotein E (APOE) is a lipid and cholesterol transport molecule known to influence Alzheimer's disease (AD) risk in an isoform-specific manner. In particular, the APOE E4 allele is the largest genetic risk factor for late-onset sporadic AD. Our recent findings uncovered activated, clonally expanded T cells in AD cerebrospinal fluid (CSF). This T cell phenotype occurred concomitantly with altered expression of APOE in CSF monocytes. Yet, whether APOE variants differentially affect peripheral immunity systems remains unknown.

**Method:** In this study, we performed targeted immune profiling using single-cell epigenetic and transcriptomic analysis of peripheral blood mononuclear cells (PBMC). We analyzed 55 age-matched healthy control (HC) and AD patients with equal distribution of APOE E3/E3, E3/E4, and E4/E4 genotypes.

**Result:** We reveal dysregulation in monocytes and clonally expanded T cells that are distinct to AD patients carrying the APOE E4/E4 genotype. Additionally, we find APOE isoform-dependent chromatin accessibility differences that correspond to RNA expression changes.

**Conclusion:** Cumulatively, these results uncover APOE isoform-dependent changes to peripheral immunity in AD.

## About
This repository contains code used to process and analyze scRNA+TCR+BCRseq and scATACseq data from the **Peripheral immunomodulatory effects of Apolipoprotein E4 in Alzheimer’s disease** study. 

All ```R``` dependencies are listed in the ```renv.lock``` file. Upon opening the ```ad_apoe_pub.Rproj``` file in ```Rstudio```, the ```renv``` package should automatically download and activate. Afterwards, the user can run ```renv::restore()``` to download all necessary packages for the study. 

The RNA+BCR+TCR dataset can be viewed and analyzed interactively in our [modified ShinyCell app](https://gatelabnu.shinyapps.io/ad_apoe_rna/).

Raw RNA+BCR+TCR .fastq files and gene expression matrices can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/). Raw ATAC .fastq files and peak matrices can be downloaded from [GSE226267](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226267).

For any comments/questions about the study please refer to the [Gate Lab website](https://sites.northwestern.edu/gatelab/) for contact information.

![Feinberg-linear-RGB](https://user-images.githubusercontent.com/91904251/221924737-8ff64f66-bc81-4155-94a3-05121b393bfc.png)


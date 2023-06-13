# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-15-2022
# Written by: Abhi Ramakrishnan, modified from Natalie Piehl's script
# Summary: Generate coverage plot for region of interest
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("GenomicRanges")
})

# set output directory
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/figures/coverageplot/out_AR_05-25-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Specify celltype of interest
celltype <- "CD8+_T_Cells"
# Degine DA range of interest 
# comment out range_oi in function if region highlight not necessary
range_oi <-  GRanges(seqnames="chr2", ranges= "127080928-127082176")
# Specify gene of interest (change to range_oi if only plotting region of interest)
gene <-  range_oi
# Specify APOE genotypes to subset
apoe <- c(
  # "E3/E3",
  # "E3/E4",
  "E4/E4"
  )
# Specify grouping variable
group <- "Diagnosis"
#-------------------------------------------------------------------------------
# Define function
custom_cov_plot <- function(celltype, gene, range_oi, apoe, group) {
  # read in seurat object
  if (celltype == "CD4+_T_Cells") {
    s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds"))
  } else {
    s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"))
  }
  
  # Map broad celltypes to object
  s@meta.data$broad_celltype <- mapvalues(s@meta.data$predictedGroupNew,
                                          from = c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                                                   "CD14_Mono", "CD16_Mono",
                                                   "ASDC", "cDC1", "cDC2", "pDC",
                                                   "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                                                   "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                                                   "NK", "NK_Proliferating", "NK_CD56bright",
                                                   "ILC", "dnT", "gdT",
                                                   "Platelet", "Eryth", "HSPC", "Doublet"),
                                          to = c(rep("B_Cells", 4),
                                                 rep("Monocytes", 2),
                                                 rep("Dendritic_Cells", 4),
                                                 rep("CD4+_T_Cells", 6),
                                                 rep("CD8+_T_Cells", 5),
                                                 rep("NK_Cells", 3),
                                                 rep("Other_T_Cells", 3),
                                                 rep("Other", 4)))
  
  # Subset for celltype of interest
  s <- subset(s, broad_celltype == celltype)
  
  # Make Sex + Diagnosis column
  s@meta.data$Diagnosis_Sex <- paste(s@meta.data$Diagnosis, s@meta.data$Sex)
  
  # Define levels
  s@meta.data$Diagnosis_Sex <- factor(s@meta.data$Diagnosis_Sex,
                                       levels = c("Healthy Control female",
                                                  "Alzheimers Disease female",
                                                  "Healthy Control male",
                                                  "Alzheimers Disease male"))
  s@meta.data$Diagnosis <- factor(s@meta.data$Diagnosis,
                                  levels = c("Healthy Control",
                                             "Alzheimers Disease"))
  s@meta.data$APOE_genotype <- factor(s@meta.data$APOE_genotype,
                                      levels = c("E3/E3",
                                                 "E3/E4",
                                                 "E4/E4"))
  
  # Downsample so each group has same number of cells
  barcodes <- s[["Diagnosis_Sex"]]
  barcodes$barcode <- rownames(barcodes)
  freq <- table(barcodes$Diagnosis_Sex) %>% as.data.frame()
  sample_num <- min(freq$Freq)
  barcodes_sampled <- barcodes %>%
    group_by(Diagnosis_Sex) %>%
    slice_sample(n = sample_num, replace = FALSE)
  s <- subset(s, cells = barcodes_sampled$barcode)
  
  # Subset for group of interest
  s <- subset(s, APOE_genotype %in% apoe)
  
  # Set identity
  Idents(s) <- group
  
  # Generate coverage plot
  cov_plot <- CoveragePlot(
    object = s,
    region = gene,
    region.highlight = range_oi,
    extend.downstream = 1000,
    extend.upstream = 1000,
    annotation = TRUE,
    peaks = TRUE,
    # ymax = 100
  )
  # redefine variables for naming
  apoe_gen <- paste(apoe, collapse = "_")
  apoe_gen <- gsub("/", "-", apoe_gen)
  range_oi <- gsub(":", "-", range_oi)
  
  pdf(paste0(output_dir, celltype, "_", gene,"_",range_oi, "_",group,"_", apoe_gen,"_covplot.pdf"))
  print(cov_plot)
  dev.off()
}

# Run function 
custom_cov_plot(celltype, gene, range_oi, apoe, group)

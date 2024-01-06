# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-10-2023
# Written by: Abhi Ramakrishnan
# Summary: Calculate frequencies of cell types across RNA-seq and ATAC-seq data 
# by diagnosis and APOE Genotype for AD-APOE project
#
# ------------------------------------------------------------------------------

# renv::restore()

# Create output directory
output_dir <- "/path/to/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# Source Natalie's helper functions
source("/path/to/helper_functions.R")

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("dplyr")
})

# Function to omit rows with zeros
omit_zero <- function(data){
  data[data==0] <- NA
  data <- data[complete.cases(data),]
}

#-------------------------------------------------------------------------------
# Part 1: All cells for fine celltypes in RNA data
#-------------------------------------------------------------------------------
# Load object with B and T cell data
s <- "/path/to/seurat/object/"
load(s)

# Extract relevant columns from matrix
data <- as.data.frame(cbind(s$orig.ident, s$predicted.celltype.l2, s$Diagnosis, s$APOE_genotype, s$Sex, s$Age, s$T_frequency_filtered))
colnames(data) <- c("Sample_ID", "celltype", "Diagnosis", "APOE_genotype", "Sex", "Age", "T_clone_size")
# Add disease_genotype column
data$DiagXApoe <- paste(data$Diagnosis, data$APOE_genotype, sep = "_")
data$DiagXApoe <- factor(data$DiagXApoe, levels = c("Healthy Control_E3/E3",
                                                    "Healthy Control_E3/E4",
                                                    "Healthy Control_E4/E4",
                                                    "Alzheimers Disease_E3/E3",
                                                    "Alzheimers Disease_E3/E4",
                                                    "Alzheimers Disease_E4/E4"))
print("Test 1")
head(data)
# #-------------------------------------------------------------------------------
# 1.) Get per sample average cell number, and total number of cells per group
# #-------------------------------------------------------------------------------
# a. where group is DiagnosisxAPOExCelltype--175 total groups
diag_ex <- data.frame(table(data$Diagnosis, data$APOE_genotype, data$celltype, data$Sample_ID))
diag_ex <- omit_zero(diag_ex)
colnames(diag_ex) <- c("Diagnosis","APOE_genotype","Celltype", "Sample_ID", "Freq")
sum_stats <- diag_ex %>%
  dplyr::group_by(Diagnosis, APOE_genotype, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                       " (",
                                       quantile(Freq, na.rm = TRUE)[2],
                                       " - ",
                                       quantile(Freq, na.rm = TRUE)[4],
                                       ")"),
                   total_cells_per_group = sum(Freq))

write.csv(sum_stats, file = paste0(output_dir, "rna_diagxGenoxCelltype_cell_counts.csv"))
# b. where group is DiagnosisxCelltype
sum_stats <- diag_ex %>%
  dplyr::group_by(Diagnosis, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))
write.csv(sum_stats, file = paste0(output_dir, "rna_diagxCelltype_cell_counts.csv"))
# #-------------------------------------------------------------------------------
# 2.) Repeat above for CLONAL cells only; Get per sample average cell number, and total number of cells per group
# #-------------------------------------------------------------------------------
# First filter for clonal cells only
data_clonal <- drop_na(data, T_clone_size)
data_clonal$T_clone_size <- as.numeric(data_clonal$T_clone_size)
head(data_clonal)
data_clonal <- dplyr::filter(data_clonal, T_clone_size > 1)
print("After filtering for clonal cells: ")
head(data_clonal)

# a. where group is DiagnosisxAPOExCelltype--175 total groups
diag_ex <- data.frame(table(data_clonal$Diagnosis, data_clonal$APOE_genotype, data_clonal$celltype, data_clonal$Sample_ID))
diag_ex <- omit_zero(diag_ex)
colnames(diag_ex) <- c("Diagnosis","APOE_genotype","Celltype", "Sample_ID", "Freq")
sum_stats <- diag_ex %>%
  dplyr::group_by(Diagnosis, APOE_genotype, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))

write.csv(sum_stats, file = paste0(output_dir, "rna_clonal-T_diagxGenoxCelltype_cell_counts.csv"))
# b. where group is DiagnosisxCelltype
sum_stats <- diag_ex %>%
  dplyr::group_by(Diagnosis, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))
write.csv(sum_stats, file = paste0(output_dir, "rna_clonal-T_diagxCelltype_cell_counts.csv"))


#-------------------------------------------------------------------------------
# Part 2: All cells for broad celltypes in ATAC data
#-------------------------------------------------------------------------------

## A.) Load CD4 object
s <- readRDS(paste0("/path/to/seurat/object"))
s[[]]
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

# Count cells

atac <- data.frame(table(s@meta.data$Diagnosis,
                         s@meta.data$APOE_genotype,
                         s@meta.data$broad_celltype,
                         s@meta.data$Sample))

# a. where group is DiagnosisxAPOExCelltype
atac <- omit_zero(atac)
colnames(atac) <- c("Diagnosis","APOE_genotype","Celltype", "Sample_ID", "Freq")
sum_stats <- atac %>%
  dplyr::group_by(Diagnosis, APOE_genotype, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                       " (",
                                       quantile(Freq, na.rm = TRUE)[2],
                                       " - ",
                                       quantile(Freq, na.rm = TRUE)[4],
                                       ")"),
                   total_cells_per_group = sum(Freq))

write.csv(sum_stats, file = paste0(output_dir, "cd4_atac_diagxGenoxCelltype_cell_counts.csv"))
# b. where group is DiagnosisxCelltype
sum_stats <- atac %>%
  dplyr::group_by(Diagnosis, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))
write.csv(sum_stats, file = paste0(output_dir, "cd4_atac_diagxCelltype_cell_counts.csv"))

## B.) Load nonCD4 object
s <- readRDS(paste0("/path/to/seurat/object"))
s[[]]
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

# Count cells

atac <- data.frame(table(s@meta.data$Diagnosis,
                         s@meta.data$APOE_genotype,
                         s@meta.data$broad_celltype,
                         s@meta.data$Sample))

# a. where group is DiagnosisxAPOExCelltype
atac <- omit_zero(atac)
colnames(atac) <- c("Diagnosis","APOE_genotype","Celltype", "Sample_ID", "Freq")
sum_stats <- atac %>%
  dplyr::group_by(Diagnosis, APOE_genotype, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))

write.csv(sum_stats, file = paste0(output_dir, "noncd4_atac_diagxGenoxCelltype_cell_counts.csv"))
# b. where group is DiagnosisxCelltype
sum_stats <- atac %>%
  dplyr::group_by(Diagnosis, Celltype) %>%
  dplyr::summarize(median_IQR_by_sample = paste0(median(Freq, na.rm = TRUE),
                                                 " (",
                                                 quantile(Freq, na.rm = TRUE)[2],
                                                 " - ",
                                                 quantile(Freq, na.rm = TRUE)[4],
                                                 ")"),
                   total_cells_per_group = sum(Freq))
write.csv(sum_stats, file = paste0(output_dir, "noncd4_atac_diagxCelltype_cell_counts.csv"))


# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-13-2022
# Written by: Morabito, Natalie Piehl
# Summary: Follow Morabito code to correlate RNA expression and ATAC
# accessibility
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("optparse")
  library("pbapply")
})

# Save arguments to variable
atac_cell_type <- "Monocytes"
rna_cell_type <- "CD14 Mono"

# Organize inputs
rna_seurat_object <-"/projects/b1169/projects/AD_APOE/results/seurat/supervised_clustering/out_NP_09-07-2022/s_sup_clustering"
cicero_dir <- "/projects/b1169/projects/AD_APOE/results_atac/cicero/batch/out_NP_02-20-2023/"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/expression_accessibility_correlation/ABCA1_exp_acc_corr/out_NP_06-08-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Prep
#-------------------------------------------------------------------------------

# Load RNA seurat object and ATAC proj + seurat object
load(rna_seurat_object)

# Subset for celltype
s <- subset(s, predicted.celltype.l2 == rna_cell_type)

if (atac_cell_type == "CD4+_T_Cells") {
  atac_s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds"))
} else {
  atac_s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"))
}

# Add broad cell types
atac_s@meta.data$broad_celltype <- mapvalues(atac_s@meta.data$predictedGroupNew,
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

# Subset for celltype
atac_s <- subset(atac_s, broad_celltype == atac_cell_type)

# Subset RNA and ATAC seurat objects by celltype and reserve only shared samples
shared_samples <- intersect(unique(s[["orig.ident"]])[,1],
                            sapply(unique(atac_s[["Sample"]])[,1],
                                   function(x) {gsub("A", "G", x)}))
s <- subset(s, orig.ident %in% shared_samples)
atac_s <- subset(atac_s, Sample %in% sapply(shared_samples,
                                            function(x) {gsub("G", "A", x)}))

# Subset for AD
s <- subset(s, Diagnosis == "Alzheimers Disease")
atac_s <- subset(atac_s, Diagnosis == "Alzheimers Disease")

# Set idents to sampleIDs:
Idents(s) <- "orig.ident"
Idents(atac_s) <- "Sample"
DefaultAssay(s) <- "RNA"
DefaultAssay(atac_s) <- "peaks"

#-------------------------------------------------------------------------------
# Define functions
#-------------------------------------------------------------------------------

# Define gene and region
gene <- "ABCA1"
region <- "chr9-104909664-104910164"

# Compute average exp/acc
avg_acc <- AverageExpression(atac_s)$peaks
avg_exp <- AverageExpression(s)$RNA

# Ensure samples are in same order
atac_sample_order <- sapply(colnames(avg_acc), function(x) {gsub("A", "G", x)})
avg_exp <- avg_exp[,atac_sample_order]

# Grab average exp and acc for this gene and peak
avg_acc_gene <- avg_acc[as.character(region),] %>% as.vector
avg_exp_gene <- avg_exp[gene,] %>% as.vector

# Grab diagnosis metadata and merge
df <- data.frame(samples = colnames(avg_acc),
                    acc = avg_acc_gene,
                    exp = avg_exp_gene)
acc <- df$acc
exp <- df$exp

# Generate ggplot object
p <- ggplot(df, aes(x = acc, y = exp)) +
  theme_Publication_blank() +
  geom_point(size = 2, shape = 21, alpha = 0.7, fill = "red") +
  # stat_cor(method = "pearson") +
  geom_smooth(method = 'lm', alpha = 0.05, linetype = 2, color = "red") +
  theme(aspect.ratio = 1,
        text = element_text(size=16))
print(p)

set_panel_size(p, file = paste0(output_dir, "ABCA1_exp_acc_correlation_scatter.pdf"),
               width = unit(4, "inch"), height = unit(4, "inch"))
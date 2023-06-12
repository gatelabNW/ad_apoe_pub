# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-18-2023
# Written by: Morabito, Natalie Piehl
# Summary: Identify gl-CREs + DARs and cre-linked genes + DEGs
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ggridges")
  library("UpSetR")
  library("stringr")
})

# Organize inputs
ranges_path <- "/projects/b1169/projects/AD_APOE/data/ranges/full_ranges.rds"
de_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis_withEthnicity/out_NP_05-18-2023_covarSex+APOE+Race//"
da_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/main/out_NP_02-06-2023/ADvsHC/"
da_intersection_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/LR_DESeq2/out_NP_05-08-2023/full_intersection/"
de_intersection_dir <- "/projects/b1169/projects/AD_APOE/results/de/MAST_edgeR/out_NP_05-19-2023/full_intersection/"
cicero_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/cicero/batch/out_NP_02-20-2023/"
input_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/expression_accessibility_correlation/batch/out_NP_03-15-2023/"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/expression_accessibility_correlation/dar_deg_comparison/out_NP_05-19-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Load in results
#-------------------------------------------------------------------------------

# Define cell type comparisons
# comparisons <- c("CD8+_T_Cells-CD8_TCM",
#                  "B_Cells-B_naive",
#                  "Monocytes-CD14_Mono",
#                  "Monocytes-CD16_Mono")
comparisons <- list.dirs(input_base_dir, full.names = FALSE)[-c(1)]

# Define diagnoses
diagnoses <- c("Healthy_Control", "Alzheimers_Disease")

# Initialize lists
glCRE_list <- c()
CRElinkedgenes_list <- c()

# Iterate through comparisons
for (comparison in comparisons) {
  # Iterate through diagnoses
  for (diagnosis in diagnoses) {
    tryCatch({
      # Read in results
      df <- readRDS(paste0(input_base_dir, comparison, "/", diagnosis, "__", comparison,
                           '_peak_gene_correlation.rds'))
      
      # Subset for 95th percentile of PCC and adjusted p val < .05
      pcc_thresh <- quantile(abs(df$pcc), 0.95)
      padj_thresh <- .05
      df <- df[which(abs(df$pcc) >= pcc_thresh & df$BH <= padj_thresh),]
      
      # Change commas to spaces
      df$avg_exp <- sapply(df$avg_exp, function(x) {gsub(",", " ", x)})
      df$avg_acc <- sapply(df$avg_acc, function(x) {gsub(",", " ", x)})
      
      # Export sig correlations
      write.table(df,
                  file=paste0(output_dir, diagnosis, "__", 
                              comparison,
                              '_peak_gene_sig_correlation.csv'),
                  sep=',', quote=FALSE, row.names=FALSE)
      
      # Isolate glCREs and cre-linked genes
      glCRE_list[[comparison]][[diagnosis]] <- unique(df$Peak2)
      CRElinkedgenes_list[[comparison]][[diagnosis]] <- unique(df$Peak1_nearestGene)
    }, error = function(e) {
      glCRE_list[[comparison]][[diagnosis]] <- c()
      CRElinkedgenes_list[[comparison]][[diagnosis]] <- c()
    })
  }
}

#-------------------------------------------------------------------------------
# Make upset of all, HC, and AD + ADvsHC degs
#-------------------------------------------------------------------------------

# Read in shared DEGs
shared_degs <- read.csv(paste0(de_intersection_dir, "ADvsHC_padj0.05_lfc0.125_MAST+edgeR_DEGs.csv"))

for (comparison in comparisons) {
  tryCatch({
    print(comparison)
    # Subset for cre-linked genes
    comp_genes <- CRElinkedgenes_list[[comparison]]
    
    # Subset for groups
    sig_genes_ls <- list()
    sig_genes_ls[["AD cre-linked genes"]] <- comp_genes[["Alzheimers_Disease"]]
    sig_genes_ls[["HC cre-linked genes"]] <- comp_genes[["Healthy_Control"]]

    # Extract rna cell type and load in degs
    rna_cell_type <- str_replace(comparison, '.+-(.+)', '\\1')
    padj.thresh <- 0.05
    lfc.thresh <- 0.125
    degs <- read.csv(paste0(de_dir, rna_cell_type, "_AD_vs_HC_degs.csv"))
    
    # Subset for shared with edgeR
    degs <- degs[which(degs$X %in% shared_degs[,gsub("\\+", ".", rna_cell_type)]),]
    
    # Select up and down degs
    up_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC > lfc.thresh),]
    down_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC < -lfc.thresh),]
    sig_genes_ls[["Up ADvsHC DEGs"]] <- up_degs$X
    sig_genes_ls[["Down ADvsHC DEGs"]] <- down_degs$X
    
    # Define color map
    color_map <- data.frame(analysis = c("HC cre-linked genes", "AD cre-linked genes", "Up ADvsHC DEGs", "Down ADvsHC DEGs"),
                            color = c("gray", "red", "turquoise1", "lawngreen"))
    
    # Order color map by number of genes
    # Find number of genes in each set
    num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
    for (key in names(sig_genes_ls)) {
      num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
    }
    
    # Get order of sets for coloring
    num_degs <- num_degs[order(-num_degs[, 2]),]
    color_map <- color_map[match(num_degs[,1], color_map$analysis),]
    
    # Create upset plot and export
    pdf(file = paste0(output_dir, "/", comparison, "_upset_ADvsHC_crelinkedGenes.pdf"))
    print(upset(
      fromList(sig_genes_ls),
      text.scale = 2,
      nsets = length(sig_genes_ls),
      order.by = "freq",
      sets.bar.color = color_map$color
    ))
    dev.off()
  }, error = function(e) {
    message(e)
  })
}

# Extract AD specific cre-linked genes + DEGs
setdiff(intersect(sig_genes_ls[["AD cre-linked genes"]], sig_genes_ls[["Up ADvsHC DEGs"]]), sig_genes_ls[["HC cre-linked genes"]])
setdiff(intersect(sig_genes_ls[["AD cre-linked genes"]], sig_genes_ls[["Down ADvsHC DEGs"]]), sig_genes_ls[["HC cre-linked genes"]])
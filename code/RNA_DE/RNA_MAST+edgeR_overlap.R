# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-02-2023
# Written by: Natalie Piehl
# Summary: Consolidate DEG lists
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("UpSetR")
  library("ComplexHeatmap")
})

# Organize inputs
de_base_dir <- "/projects/b1169/projects/AD_APOE/results/"
output_base_dir <- "/projects/b1169/projects/AD_APOE/results/de/MAST_edgeR/out_NP_05-19-2023/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparison of interests
comparison <- "ADvsHC_44"

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Define threshold specific output dir
output_dir <- paste0(output_base_dir, comparison, "/", "padj", as.character(padj.thresh),
                     "_lfc", as.character(lfc.thresh), "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_base_dir, "full_intersection"), showWarnings = FALSE)

# #-------------------------------------------------------------------------------
# # Identify cell types to remove
# #-------------------------------------------------------------------------------
# 
# rna_seurat_object <-"/projects/b1169/projects/AD_APOE/results/seurat/supervised_clustering/out_NP_09-07-2022/s_sup_clustering"
# load(rna_seurat_object)
# 
# s@meta.data$E4_carrier <- mapvalues(s@meta.data$APOE_genotype,
#                                     from = c("E3/E3", "E3/E4", "E4/E4"),
#                                     to = c("noncarrier", "carrier", "carrier"))
# 
# s@meta.data$Diagnosis_APOE <- paste(s@meta.data$Diagnosis, s@meta.data$APOE_genotype)
# s@meta.data$Diagnosis_APOE4carrier <- paste(s@meta.data$Diagnosis, s@meta.data$E4_carrier)
# 
# diagnosis_freqs <- table(s[[c("predicted.celltype.l2", "Diagnosis")]]) %>% data.frame()
# apoe_freqs <- table(s[[c("predicted.celltype.l2", "Diagnosis_APOE")]]) %>% data.frame()
# carrier_freqs <- table(s[[c("predicted.celltype.l2", "Diagnosis_APOE4carrier")]]) %>% data.frame()
# just_carrier_freqs <- table(s[[c("predicted.celltype.l2", "E4_carrier")]]) %>% data.frame()
# 
# if (comparison == "ADvsHC") {
#   group_1 <- "Healthy Control"
#   group_2 <- "Alzheimers Disease"
# 
#   group_1_celltypes_rm <- diagnosis_freqs[which(diagnosis_freqs$Diagnosis == group_1 & diagnosis_freqs$Freq < 100), "predicted.celltype.l2"]
#   group_2_celltypes_rm <- diagnosis_freqs[which(diagnosis_freqs$Diagnosis == group_2 & diagnosis_freqs$Freq < 100), "predicted.celltype.l2"]
# } else if (comparison == "44vs33_") {
#   # Specify groups
#   group_1 <- "Healthy Control E4/E4"
#   group_2 <- "Healthy Control E3/E4"
# 
#   # Find celltypes with less than 100 cells in each group
#   group_1_celltypes_rm <- apoe_freqs[which(apoe_freqs$Diagnosis_APOE == group_1 & apoe_freqs$Freq < 100), "predicted.celltype.l2"]
#   group_2_celltypes_rm <- apoe_freqs[which(apoe_freqs$Diagnosis_APOE == group_2 & apoe_freqs$Freq < 100), "predicted.celltype.l2"]
# } else if (comparison == "4carriervs33_AD") {
#   # Specify groups
#   group_1 <- "Alzheimers Disease carrier"
#   group_2 <- "Healthy Control carrier"
# 
#   # Find celltypes with less than 100 cells in each group
#   group_1_celltypes_rm <- carrier_freqs[which(carrier_freqs$Diagnosis_APOE4carrier == group_1 & carrier_freqs$Freq < 100), "predicted.celltype.l2"]
#   group_2_celltypes_rm <- carrier_freqs[which(carrier_freqs$Diagnosis_APOE4carrier == group_2 & carrier_freqs$Freq < 100), "predicted.celltype.l2"]
#   # group_1_celltypes_rm <- just_carrier_freqs[which(just_carrier_freqs$E4_carrier == group_1 & just_carrier_freqs$Freq < 100), "predicted.celltype.l2"]
#   # group_2_celltypes_rm <- just_carrier_freqs[which(just_carrier_freqs$v == group_2 & just_carrier_freqs$Freq < 100), "predicted.celltype.l2"]
# }
# celltypes_rm <- cat(union(group_1_celltypes_rm, group_2_celltypes_rm), sep = "', '")

#-------------------------------------------------------------------------------
# Prep
#-------------------------------------------------------------------------------

# Define cell types
cell_types <- c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                "CD14_Mono", "CD16_Mono",
                "ASDC", "cDC1", "cDC2", "pDC",
                "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                "NK", "NK_Proliferating", "NK_CD56bright",
                "ILC", "dnT", "gdT")

# Define comparison deg dirs
if (comparison == "ADvsHC") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'HSPC', 'NK_Proliferating', 'pDC', 'Eryth', 'Plasmablast')
} else if (comparison == "ADvsHC_33") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} else if (comparison == "ADvsHC_44") {
  mast_dir <- paste0(de_base_dir, "de/diagnosis_44/out_NP_09-28-2022_covarSex/")
  celltypes_rm <- c('ASDC', 'CD16_Mono', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "ADvsHC_34") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs33_AD") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'CD16_Mono', 'cDC2')
} else if (comparison == "34vs33_AD") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs34_AD") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'CD16_Mono', 'Platelet')
} else if (comparison == "44vs33_HC") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} else if (comparison == "34vs33_HC") {
  mast_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs34_HC") {
  mast_dir <- paste0(de_base_dir, "de/apoe44vs34_hc/out_NP_11-10-2022_covarSex/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} 
# else if (comparison == "ADvsHC_4carrier") {
#   mast_dir <- paste0(de_base_dir, "diagnosis_4carrier/out_NP_05-01-2023_covarSex/")
#   celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'Doublet', 'Eryth', 'HSPC', 'pDC', 'Plasmablast', 'dnT', 'ILC', 'NK_Proliferating', 'Platelet')
# } else if (comparison == "4carriervs33_AD") {
#   mast_dir <- paste0(de_base_dir, "apoe4carriervs33_ad_withEthnicity/out_NP_05-18-2023_covarSex+Race/")
#   celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'Doublet', 'Eryth', 'HSPC', 'pDC', 'Plasmablast', 'dnT', 'ILC', 'NK_Proliferating', 'Platelet')
# } else if (comparison == "4carriervs33_HC") {
#   mast_dir <- paste0(de_base_dir, "apoe4carriervs33_hc_withEthnicity/out_NP_05-18-2023_covarSex+Race/")
#   celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
# } 
# else if (comparison == "4carriervs33") {
#   mast_dir <- paste0(de_base_dir, "4carriervs33/out_NP_05-05-2023_covarSexDiagnosis/")
#   celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'Doublet', 'Eryth', 'HSPC')
# } else if (comparison == "44vs33") {
#   mast_dir <- paste0(de_base_dir, "44vs33/out_NP_05-05-2023_covarSexDiagnosis/")
#   celltypes_rm <- c('')
# } else if (comparison == "34vs33") {
#   mast_dir <- paste0(de_base_dir, "34vs33/out_NP_05-05-2023_covarSexDiagnosis/")
#   celltypes_rm <- c('')
# }
pseudobulk_dir <- paste0(de_base_dir, "de/DElegate/out_NP_05-05-2023/", comparison, "/")

# Initialize deg lists
mast_up_list <- list()
mast_down_list <- list()
pseudobulk_up_list <- list()
pseudobulk_down_list <- list()
down_list <- list()
up_list <- list()
combined_list <- list()

# Define deg function
grab_degs <- function(cell_type, de_dir, direction) {
  tryCatch({
    # Read in DEGs
    celltype_files <- list.files(de_dir, pattern = cell_type, full.names = TRUE)
    if (cell_type == "NK") {
      files <- grep(".csv", celltype_files, value = TRUE)
      degs <- read.csv(files[!grepl("CD56|Prolif", files)])
    } else {
      degs <- read.csv(grep(".csv", celltype_files, value = TRUE))
    }
    
    if (de_dir == pseudobulk_dir) {
      # Identify sig genes
      up_degs <- degs[which(degs$pvalue < padj.thresh & degs$log_fc > lfc.thresh),]
      down_degs <- degs[which(degs$pvalue < padj.thresh & degs$log_fc < -lfc.thresh),]
      
      # Save genes to list
      if (direction == "up") {
        return(up_degs$feature)
      } else {return(down_degs$feature)}
    } else {
      # Identify sig genes
      up_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC > lfc.thresh),]
      down_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC < -lfc.thresh),]
      
      # Save genes to list
      if (direction == "up") {
        return(up_degs$X)
      } else {return(down_degs$X)}
    }
  }, error = function(e) {
    return(c())
  })
}

# Initialize MAST and pseudobulk overlap list
mast_pseudobulk_list <- list()

#-------------------------------------------------------------------------------
# Generate upsets per cell type
#-------------------------------------------------------------------------------

for (cell_type in cell_types[ cell_types %!in% celltypes_rm ]) {
  # Grab MAST DEGs
  mast_up_list[[cell_type]] <- grab_degs(cell_type, mast_dir, "up")
  mast_down_list[[cell_type]] <- grab_degs(cell_type, mast_dir, "down")
  
  # Grab edgeR DEGs
  pseudobulk_up_list[[cell_type]] <- grab_degs(cell_type, pseudobulk_dir, "up")
  pseudobulk_down_list[[cell_type]] <- grab_degs(cell_type, pseudobulk_dir, "down")

  # Generate up and down lists
  down_list[[cell_type]] <- list("MAST" = mast_down_list[[cell_type]],
                                 "Pseudobulk" = pseudobulk_down_list[[cell_type]])
  up_list[[cell_type]] <- list("MAST" = mast_up_list[[cell_type]],
                                 "Pseudobulk" = pseudobulk_up_list[[cell_type]])

  # Define color map
  color_map <- data.frame(method = c("MAST", "Pseudobulk"),
                          color = c("cyan3", "plum2"))

  # Retain only nonzero groups
  down_list[[cell_type]] <- down_list[[cell_type]][ lapply(down_list[[cell_type]], length) > 0 ]
  up_list[[cell_type]] <- up_list[[cell_type]][ lapply(up_list[[cell_type]], length) > 0 ]

  # Up DEGs --------------------------------------------------------------------
  if (length(up_list[[cell_type]]) > 1) {
    tryCatch({
      # Save results
      intersected_degs <- intersect(up_list[[cell_type]][["MAST"]], up_list[[cell_type]][["Pseudobulk"]])
      write.csv(intersected_degs, paste0(output_dir, comparison, "_", cell_type, "_Up-MAST+edgeR_DEGs.csv"))
    }, error = function(e) {
      message(e)
    })
  }

  # Down DEGs ------------------------------------------------------------------
  if (length(down_list[[cell_type]]) > 1) {
    tryCatch({
      # Save results
      intersected_degs <- intersect(down_list[[cell_type]][["MAST"]], down_list[[cell_type]][["Pseudobulk"]])
      write.csv(intersected_degs, paste0(output_dir, comparison, "_", cell_type, "_Down-MAST+edgeR_DEGs.csv"))
    }, error = function(e) {
      message(e)
    })
  }

  # Combined DEGs --------------------------------------------------------------
  # Check for discordant genes
  discordant_genes <- union(intersect(down_list[[cell_type]][["MAST"]], up_list[[cell_type]][["MAST"]]),
                            intersect(down_list[[cell_type]][["Pseudobulk"]], up_list[[cell_type]][["Pseudobulk"]]))
  if (length(discordant_genes) != 0) {
    print(discordant_genes)
    next
  }

  # Make combined DEG lists
  combined_list[[cell_type]] <- list("MAST" = union(down_list[[cell_type]][["MAST"]], up_list[[cell_type]][["MAST"]]),
                                     "Pseudobulk" = union(down_list[[cell_type]][["Pseudobulk"]], up_list[[cell_type]][["Pseudobulk"]]))

  # Retain only nonzero groups
  combined_list[[cell_type]] <- combined_list[[cell_type]][ lapply(combined_list[[cell_type]], length) > 0 ]

  if (length(combined_list[[cell_type]]) > 1) {
    tryCatch({
      # Find number of genes in each set
      num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
      for (key in names(combined_list[[cell_type]])) {
        num_degs <- rbind(num_degs, data.frame(key, length(combined_list[[cell_type]][[key]])))
      }

      # Get order of sets for coloring
      num_degs <- num_degs[order(-num_degs[, 2]),]
      color_map <- color_map[match(num_degs[,1], color_map$method),]

      # Create upset plot and export
      pdf(file = paste0(output_dir, comparison, "_", cell_type, "_DEGoverlap_upset.pdf"))
      print(upset(
        fromList(combined_list[[cell_type]]),
        text.scale = 3,
        nsets = length(combined_list[[cell_type]]),
        order.by = "freq",
        sets.bar.color = color_map$color
      ))
      dev.off()

      # Add intersection to list
      mast_pseudobulk_list[[cell_type]] <- intersect(combined_list[[cell_type]][["MAST"]],
                                                     combined_list[[cell_type]][["Pseudobulk"]])
    }, error = function(e) {
      message(e)
    })
  }
}

# Export intersected list
max_degs <- 0
for (cell_type in cell_types[ cell_types %!in% celltypes_rm ]) {
  degs <- mast_pseudobulk_list[[cell_type]]
  if (length(degs) > max_degs) {
    max_degs <- length(degs) 
  }
}
col_list <- list()
for (cell_type in cell_types[ cell_types %!in% celltypes_rm ]) {
  degs <- mast_pseudobulk_list[[cell_type]]
  col_list[[cell_type]] <- c(degs, rep(NA, max_degs - length(degs)))
}
intersected_df <- do.call(cbind, col_list) %>% data.frame
intersected_df[is.na(intersected_df)] <- ""
write.csv(intersected_df, paste0(output_base_dir, "full_intersection/", comparison, "_",
                                 "padj", as.character(padj.thresh),
                                 "_lfc", as.character(lfc.thresh),
                                 "_MAST+edgeR_DEGs.csv"))
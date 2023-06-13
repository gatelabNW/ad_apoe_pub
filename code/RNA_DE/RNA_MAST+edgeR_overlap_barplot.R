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
  library("UpSetR")
  library("ComplexHeatmap")
})

# Organize inputs
celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/celltype_color_map.csv"
de_base_dir <- "/projects/b1169/projects/AD_APOE/results/"
output_dir <- "/projects/b1169/projects/AD_APOE/results/de/MAST_edgeR_allcelltypes_upset/out_NP_05-24-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparison of interests
comparison <- "ADvsHC"

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

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
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'HSPC', 'NK_Proliferating', 'pDC', 'Eryth', 'Plasmablast')
} else if (comparison == "ADvsHC_33") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} else if (comparison == "ADvsHC_44") {
  MAST_dir <- paste0(de_base_dir, "de/diagnosis_44/out_NP_09-28-2022_covarSex/")
  celltypes_rm <- c('ASDC', 'CD16_Mono', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "ADvsHC_34") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs33_AD") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'CD16_Mono', 'cDC2')
} else if (comparison == "34vs33_AD") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs34_AD") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'CD16_Mono', 'Platelet')
} else if (comparison == "44vs33_HC") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} else if (comparison == "34vs33_HC") {
  MAST_dir <- paste0(de_base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
} else if (comparison == "44vs34_HC") {
  MAST_dir <- paste0(de_base_dir, "de/apoe44vs34_hc/out_NP_11-10-2022_covarSex/")
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
} 
pseudobulk_dir <- paste0(de_base_dir, "de/DElegate/out_NP_05-05-2023/", comparison, "/")

# Initialize deg lists
MAST_up_list <- list()
MAST_down_list <- list()
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
MAST_pseudobulk_list <- list()

#-------------------------------------------------------------------------------
# Generate upsets per cell type
#-------------------------------------------------------------------------------

grab_overlaps <- function(cell_type) {
  # Grab MAST DEGs
  MAST_up_list[[cell_type]] <- grab_degs(cell_type, MAST_dir, "up")
  MAST_down_list[[cell_type]] <- grab_degs(cell_type, MAST_dir, "down")
  MAST_genes <- union(MAST_up_list[[cell_type]], MAST_down_list[[cell_type]])
  
  # Grab edgeR DEGs
  pseudobulk_up_list[[cell_type]] <- grab_degs(cell_type, pseudobulk_dir, "up")
  pseudobulk_down_list[[cell_type]] <- grab_degs(cell_type, pseudobulk_dir, "down")
  pseudobulk_genes <- union(pseudobulk_up_list[[cell_type]], pseudobulk_down_list[[cell_type]])

  # Grab num unique and shared
  MAST_unique <- MAST_genes %>% length
  pb_unique <- pseudobulk_genes %>% length
  shared <- intersect(pseudobulk_genes, MAST_genes) %>% length
  
  return(c(MAST_unique, pb_unique, shared))
}

# Calculate overlaps
res <- lapply(cell_types[cell_types %!in% celltypes_rm], grab_overlaps)
df <- do.call(rbind, res) %>% as.data.frame
colnames(df) <- c("MAST only", "edgeR only", "shared")
rownames(df) <- cell_types[cell_types %!in% celltypes_rm]

# Convert to long format
df$cell_type <- rownames(df)
df_long <- df %>%
  pivot_longer(cols = -c(cell_type),
               names_to = "Comparison",
               values_to = 'Num_Degs')
cell_type_order <- rownames(df[order(df$shared, decreasing = TRUE),])
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

# Grab colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[nrow(celltype_colors) + 1,] <- c("shared", "gray", "gray")
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
# celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub("\\+", ".", x)})
celltype_colors <- celltype_colors[match(cell_type_order, celltype_colors$predicted.celltype.l2),]
df_long$Comparison <- factor(df_long$Comparison, levels = c("shared", "edgeR only", "MAST only"))

# Generate plot
p <- ggplot(df_long, aes(y = Comparison, x = Num_Degs, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = celltype_colors$new_color) +
  theme_Publication_blank() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right",
        text = element_text(size = 24))
p

set_panel_size(p, file = paste0(output_dir, "ADvsHC_MASTvsedgeR_TotalnumDEGs_celltype_origin.pdf"),
               width = unit(7, "in"), height = unit(4, "in"))
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
# Summary: Consolidate DAR lists
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
da_celltype_colors_path <- "/path/to/broad_celltype_color_map.csv"
da_base_dir <- "path/to/DA/results/"
celltype_colors_path <- "/path/to/fine_celltype_color_map.csv"
output_dir <- "/path/to/output/dir/"
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
cell_types <- c("B_Cells", "Monocytes", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells")

# Define comparison deg dirs
LR_dir <- paste0(da_base_dir, "main/out_NP_02-06-2023/", comparison, "/")
pseudobulk_dir <- paste0(da_base_dir, "DElegate/out_NP_05-05-2023/", comparison, "/")

# Initialize deg lists
LR_up_list <- list()
LR_down_list <- list()
pseudobulk_up_list <- list()
pseudobulk_down_list <- list()
down_list <- list()
up_list <- list()
combined_list <- list()

# Define deg function
grab_dars <- function(cell_type, de_dir, direction) {
  tryCatch({
    # Read in DEGs
    if (cell_type == "CD8+_T_Cells") {
      celltype_files <- list.files(de_dir, pattern = "CD8", full.names = TRUE)
    } else if (cell_type == "CD4+_T_Cells") {
      celltype_files <- list.files(de_dir, pattern = "CD4", full.names = TRUE)
    } else {
      celltype_files <- list.files(de_dir, pattern = cell_type, full.names = TRUE)
    }
    degs <- read.csv(grep(".csv", celltype_files, value = TRUE))

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
        return(up_degs$sitename)
      } else {return(down_degs$sitename)}
    }
  }, error = function(e) {
    return(c())
  })
}

# Initialize LR and pseudobulk overlap list
LR_pseudobulk_list <- list()

#-------------------------------------------------------------------------------
# Generate upsets per cell type
#-------------------------------------------------------------------------------

grab_overlaps <- function(cell_type) {
  # Grab LR DEGs
  LR_up_list[[cell_type]] <- grab_dars(cell_type, LR_dir, "up")
  LR_down_list[[cell_type]] <- grab_dars(cell_type, LR_dir, "down")
  LR_genes <- union(LR_up_list[[cell_type]], LR_down_list[[cell_type]])
  
  # Grab DESeq2 DEGs
  pseudobulk_up_list[[cell_type]] <- grab_dars(cell_type, pseudobulk_dir, "up")
  pseudobulk_down_list[[cell_type]] <- grab_dars(cell_type, pseudobulk_dir, "down")
  pseudobulk_genes <- union(pseudobulk_up_list[[cell_type]], pseudobulk_down_list[[cell_type]])

  # Grab num unique and shared
  lr_unique <- LR_genes %>% length
  pb_unique <- pseudobulk_genes %>% length
  shared <- intersect(pseudobulk_genes, LR_genes) %>% length
  
  return(c(lr_unique, pb_unique, shared))
}

# Calculate overlaps
res <- lapply(cell_types, grab_overlaps)
df <- do.call(rbind, res) %>% as.data.frame
colnames(df) <- c("LR only", "DESeq2 only", "shared")
rownames(df) <- cell_types

# Convert to long format
df$cell_type <- rownames(df)
df_long <- df %>%
  pivot_longer(cols = -c(cell_type),
               names_to = "Comparison",
               values_to = 'Num_Dars')
cell_type_order <- rownames(df[order(df$shared, decreasing = TRUE),])
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

# Grab colors
celltype_colors <- read.csv(da_celltype_colors_path)
celltype_colors[nrow(celltype_colors) + 1,] <- c("shared", "gray")
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
celltype_colors <- celltype_colors[match(cell_type_order, celltype_colors$predicted.celltype.l2),]
df_long$Comparison <- factor(df_long$Comparison, levels = c("shared", "DESeq2 only", "LR only"))

# Generate plot
p <- ggplot(df_long, aes(y = Comparison, x = Num_Dars, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = celltype_colors$new_color) +
  theme_Publication_blank() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right",
        text = element_text(size = 24))
p

set_panel_size(p, file = paste0(output_dir, "ADvsHC_LRvsDESeq2_TotalnumDARs_celltype_origin.pdf"),
               width = unit(7, "in"), height = unit(4, "in"))

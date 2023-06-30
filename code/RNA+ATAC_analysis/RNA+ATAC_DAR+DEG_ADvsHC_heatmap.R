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
celltype_colors_path <- "/path/to/celltype_color_map.csv"
broad_celltype_colors_path <- "/path/to/broad_celltype_color_map.csv"
ranges_path <- "/path/to/ranges/object"
base_dir <- "/path/to/project/root/"
output_dir <- "/path/to/output_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparison of interests
comparison <- "ADvsHC"

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

#-------------------------------------------------------------------------------
# Prep
#-------------------------------------------------------------------------------

# Define celltypes to remove
if (comparison == "ADvsHC") {
  celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'HSPC', 'NK_Proliferating', 'pDC', 'Eryth', 'Plasmablast')
}

# Make broad cell type map
broad_map <- data.frame(fine = c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                                 "CD14_Mono", "CD16_Mono",
                                 "ASDC", "cDC1", "cDC2", "pDC",
                                 "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                                 "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                                 "NK", "NK_Proliferating", "NK_CD56bright",
                                 "ILC", "dnT", "gdT",
                                 "Platelet", "Eryth", "HSPC", "Doublet"),
                        broad = c(rep("B_Cells", 4),
                                  rep("Monocytes", 2),
                                  rep("Dendritic_Cells", 4),
                                  rep("CD4+_T_Cells", 6),
                                  rep("CD8+_T_Cells", 5),
                                  rep("NK_Cells", 3),
                                  rep("Other_T_Cells", 3),
                                  rep("Other", 4)))

# Define cell types
cell_types <- c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                "CD14_Mono", "CD16_Mono",
                "ASDC", "cDC1", "cDC2", "pDC",
                "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                "NK", "NK_Proliferating", "NK_CD56bright",
                "ILC", "dnT", "gdT")

# Read in intersections
full_deg_list <- read.csv(paste0(base_dir, "path/to/MAST+edgeR/results", comparison, "_padj0.05_lfc0.125_MAST+edgeR_DEGs.csv"))
full_dar_list <- read.csv(paste0(base_dir, "path/to/LR+DESeq2/results", comparison, "_padj0.05_lfc0.125_LR+DESeq2_DARs.csv"))

# Read in ranges
ranges <- readRDS(ranges_path)

#-------------------------------------------------------------------------------
# Generate barplot
#-------------------------------------------------------------------------------

# Save overlap genes for upset
sig_genes_ls <- list()

grab_overlaps <- function(cell_type) {
  tryCatch({
    # Find broad celltype
    broad_cell_type <- broad_map[which(broad_map$fine == cell_type), "broad"]
    
    # Grab intersections
    deg_list <- full_deg_list[,cell_type][!is.na(full_deg_list[,cell_type]) & full_deg_list[,cell_type] != ""]
    dar_list <- full_dar_list[,gsub("\\+", ".", broad_cell_type)][!is.na(full_dar_list[,gsub("\\+", ".", broad_cell_type)]) & full_dar_list[,gsub("\\+", ".", broad_cell_type)] != ""]

    # Merge with ranges
    dar_ranges <- ranges[which(ranges$sitename %in% dar_list),]
    
    # Grab num unique and shared
    DEGs <- deg_list %>% length
    DAR_nearestGenes <- unique(dar_ranges$nearestGene) %>% length
    shared <- dar_ranges$nearestGene[which(dar_ranges$nearestGene %in% deg_list)] %>% length
    
    # Save intersected genes for upset
    # sig_genes_ls[[cell_type]] <<- intersect(deg_list, unique(dar_ranges$nearestGene))
    sig_genes_ls[[cell_type]] <<- dar_ranges$nearestGene[which(dar_ranges$nearestGene %in% deg_list)]
    
    return(c(DEGs, DAR_nearestGenes, shared))
  }, error = function(e) {
    return(c(NA, NA, NA))
  })
}

# Calculate overlaps
res <- lapply(cell_types[cell_types %!in% celltypes_rm], grab_overlaps)
df <- do.call(rbind, res) %>% as.data.frame
colnames(df) <- c("DEGs", "DARs", "shared")
rownames(df) <- cell_types[cell_types %!in% celltypes_rm]

# Sum fine cell types within a broad cell types
df$`shared sum` <- rep(0, nrow(df))
df$`DEGs sum` <- rep(0, nrow(df))
df$cell_type <- rep(NA, nrow(df))
for (broad_type in unique(broad_map$broad)) {
  # Sum num shared and num DEGs from each fine cell type
  fine_types <- broad_map[which(broad_map$broad == broad_type), "fine"]
  df[fine_types, "shared sum"] <- sum(df[fine_types, "shared"][!is.na(df[fine_types, "shared"])])
  df[fine_types, "DEGs sum"] <- sum(df[fine_types, "DEGs"][!is.na(df[fine_types, "DEGs"])])
  df[fine_types, "cell_type"] <- broad_type
}
df <- unique(df[,c(2,4,5,6)])
df <- df[!is.na(df$DARs),]

# Convert to long format
df_long <- df %>%
  pivot_longer(cols = -c(cell_type),
               names_to = "Comparison",
               values_to = 'Num')
cell_type_order <- df[order(df$`shared sum`, decreasing = TRUE), "cell_type"]
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

# Grab colors
celltype_colors <- read.csv(broad_celltype_colors_path)
celltype_colors <- celltype_colors[match(cell_type_order, celltype_colors$predicted.celltype.l2),]
df_long$Comparison <- factor(df_long$Comparison, levels = c("shared sum", "DEGs sum", "DARs"))

# Generate plot
p <- ggplot(df_long, aes(y = Comparison, x = Num, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = celltype_colors$new_color) +
  theme_Publication_blank() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right",
        text = element_text(size = 24))
p

set_panel_size(p, file = paste0(output_dir, "ADvsHC_DEGvsDAR_overlap_celltype_origin.pdf"),
               width = unit(7, "in"), height = unit(4, "in"))

#-------------------------------------------------------------------------------
# Generate upset
#-------------------------------------------------------------------------------

# Reserve cell types with at least 1 element
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
tmp_colors <- read.csv(celltype_colors_path)
tmp_colors$predicted.celltype.l2 <- sapply(tmp_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
num_degs <- num_degs[order(-num_degs[, 2]),]
tmp_colors <- tmp_colors[match(num_degs[,1], tmp_colors$predicted.celltype.l2),]
tmp_colors <- tmp_colors[ which(tmp_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "ADvsHC_DEGvsDAR_overlap_upset.pdf"))
print(upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = tmp_colors$new_color,
  text.scale = 1.25
))
dev.off()

#-------------------------------------------------------------------------------
# Generate heatmap
#-------------------------------------------------------------------------------

# Iterate through genotypes
comparisons <- c("ADvsHC_33", "ADvsHC_34", "ADvsHC_44")
apoe_comparisons <- function(comparison) {
  if (comparison == "ADvsHC_33") {
    celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet', 'cDC2')
  } else if (comparison == "ADvsHC_44") {
    celltypes_rm <- c('ASDC', 'CD16_Mono', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
  } else if (comparison == "ADvsHC_34") {
    celltypes_rm <- c('ASDC', 'CD4_Proliferating', 'CD8_Proliferating', 'cDC2', 'dnT', 'Doublet', 'Eryth', 'HSPC', 'ILC', 'NK_Proliferating', 'pDC', 'Plasmablast', 'Platelet')
  }
  
  # Read in intersections
  full_deg_list <<- read.csv(paste0(base_dir, "path/to/MAST+edgeR/results/", comparison, "_padj0.05_lfc0.125_MAST+edgeR_DEGs.csv"))
  full_dar_list <<- read.csv(paste0(base_dir, "path/to/LR+DESeq2/results/", comparison, "_padj0.05_lfc0.125_LR+DESeq2_DARs.csv"))
  
  # Calculate overlaps
  res <- lapply(cell_types, grab_overlaps)
  df <- do.call(rbind, res) %>% as.data.frame
  colnames(df) <- c("DEGs", "DAR nearest genes", "shared")
  df$cell_type <- cell_types
  df[which(df$celltype %in% celltypes_rm), "shared"] <- NA
  df$comparison <- rep(comparison, nrow(df))
  
  return(df)
}

# Calculate results
res <- lapply(comparisons, apoe_comparisons)
df <- do.call(rbind, res) %>% as.data.frame

# Widen
df_wide <- pivot_wider(df[,3:5], values_from = "shared", names_from = "comparison")
df_wide <- df_wide[-which(rowSums(df_wide[,2:4]) == 0),]
df_wide <- df_wide[rowSums(is.na(df_wide[,2:4])) != ncol(df_wide[,2:4]), ]

# Generate with hclust
pdf(paste0(output_dir, "ADvsHC_inAPOE_DAR+DEG_overlapNum_heatmap.pdf"),
    width = 5, height = 3)
data_heatmap <- t(df_wide[,2:4])
colnames(data_heatmap) <- df_wide$cell_type
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         border_color = NA
)
dev.off()

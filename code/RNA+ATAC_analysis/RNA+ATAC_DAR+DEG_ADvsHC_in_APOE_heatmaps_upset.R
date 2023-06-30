# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 11-29-2022
# Written by: Natalie Piehl
# Summary: Generate heatmaps comparing Diagnosis DE on different APOE genotypes
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("pheatmap")
})

# Organize inputs
celltype_colors_path <- "/path/to/fine_celltype_color_map.csv"
da_celltype_colors_path <- "/path/to/broad_celltype_color_map.csv"
input_dir <- "/path/to/MAST+edgeR/results/"
da_input_dir <- "/path/to/LR+DESeq2/results/"
output_dir <- "/path/to/output_dir/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

#-------------------------------------------------------------------------------
# ADvsHC in APOE # DEGs heatmap
#-------------------------------------------------------------------------------

# Define folders to compare
comparisons <- c("ADvsHC_33", "ADvsHC_34", "ADvsHC_44")
deg_list <- list()
gene_names_list <- list()
cd14_mono_list <- list()

# Define cell types
celltype_list <- c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                "CD14_Mono", "CD16_Mono",
                "ASDC", "cDC1", "cDC2", "pDC",
                "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                "NK", "NK_Proliferating", "NK_CD56bright",
                "ILC", "dnT", "gdT")

for (comparison in comparisons) {
  # Read in intersection
  intersected_df <- read.csv(paste0(input_dir, comparison, "_",
                                    "padj", as.character(padj.thresh),
                                    "_lfc", as.character(lfc.thresh),
                                    "_MAST+edgeR_DEGs.csv"), row.names = 1)
  
  # Create list
  intersected_list <- list()
  for(cell_type in names(intersected_df)) {
    intersected_list[[cell_type]] <- intersected_df[ , cell_type][ (intersected_df[ , cell_type] != "") & (!is.na(intersected_df[ , cell_type])) ]
  }

  # Initialize deg vector
  deg_vector <- c()
  gene_names <- c()

  for (cell_type in celltype_list) {
    # Add sig genes to list
    if (cell_type %in% names(intersected_list)) {
      deg_vector <- c(deg_vector, length(intersected_list[[cell_type]]))
      gene_names <- union(gene_names, intersected_list[[cell_type]])
    } else {
      deg_vector <- c(deg_vector, NA)
    }
  }
  cd14_mono_list[[comparison]] <- intersected_list[["CD14_Mono"]]

  # Add genes to list
  deg_list[[comparison]] <- deg_vector
  gene_names_list[[comparison]] <- gene_names
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- c("E3/E3", "E3/E4", "E4/E4")
data

# Remove zero rows and all NA rows
data <- data[-which(rowSums(data) == 0),]
data <- data[rowSums(is.na(data)) != ncol(data), ]

# Generate with hclust
pdf(paste0(output_dir, "ADvsHC_inAPOE_numDEGs_heatmap.pdf"),
    width = 5, height = 3)
data_heatmap <- t(data[,1:3])
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         border_color = NA
)
dev.off()

#-------------------------------------------------------------------------------
# ADvsHC in APOE Unique DEGs bar plot
#-------------------------------------------------------------------------------

# Make upset for confirmation
pdf(paste0(output_dir, "ADvsHC_inAPOE_numDEGs_upset.pdf"))
print(upset(
  fromList(gene_names_list),
  text.scale = 2,
  sets = names(gene_names_list),
  nsets = length(gene_names_list),
  sets.bar.color = c("gold", "red", "orange")
))
dev.off()

# Make cd14 mono upset for confirmation
pdf(paste0(output_dir, "CD14_Mono_ADvsHC_inAPOE_numDEGs_upset.pdf"))
print(upset(
  fromList(cd14_mono_list),
  text.scale = 2,
  sets = names(cd14_mono_list),
  nsets = length(cd14_mono_list),
  sets.bar.color = c("red", "yellow", "orange")
))
dev.off()

# Identify unique APOE DEGs
unique_gene_names_list <- list()
unique_gene_names_list[["ADvsHC_33"]] <- setdiff(gene_names_list[["ADvsHC_33"]],
                                                 union(gene_names_list[["ADvsHC_34"]], gene_names_list[["ADvsHC_44"]]))
unique_gene_names_list[["ADvsHC_34"]] <- setdiff(gene_names_list[["ADvsHC_34"]],
                                                 union(gene_names_list[["ADvsHC_33"]], gene_names_list[["ADvsHC_44"]]))
unique_gene_names_list[["ADvsHC_44"]] <- setdiff(gene_names_list[["ADvsHC_44"]],
                                                 union(gene_names_list[["ADvsHC_34"]], gene_names_list[["ADvsHC_33"]]))

# Initialize cell type counter
df <- data.frame(row.names = c(celltype_list, "shared"),
                 "ADvsHC_33" = rep(0, length(celltype_list) + 1),
                 "ADvsHC_34" = rep(0, length(celltype_list) + 1),
                 "ADvsHC_44" = rep(0, length(celltype_list) + 1))

for (comparison in comparisons) {
  # Read in intersection
  intersected_df <- read.csv(paste0(input_dir, comparison, "_",
                                    "padj", as.character(padj.thresh),
                                    "_lfc", as.character(lfc.thresh),
                                    "_MAST+edgeR_DEGs.csv"), row.names = 1)
  
  # Create list
  intersected_list <- list()
  for(cell_type in names(intersected_df)) {
    intersected_list[[cell_type]] <- intersected_df[ , cell_type][ (intersected_df[ , cell_type] != "") & (!is.na(intersected_df[ , cell_type])) ]
  }
  
  # Specify DEGs
  gene_names <- unique_gene_names_list[[comparison]]
  for (gene in gene_names) {
    # initialize result vector
    gene_celltypes <- c()
    for (cell_type in celltype_list[celltype_list %in% names(intersected_list)]) {
      if (gene %in% intersected_list[[cell_type]]) {
        gene_celltypes <- c(gene_celltypes, cell_type)
      }
    }
    if (length(gene_celltypes) == 1) {
      df[gene_celltypes[1], comparison] <- df[gene_celltypes[1], comparison] + 1
    } else {
      df["shared", comparison] <- df["shared", comparison] + 1
    }
  }
}

# Convert to long format
df$cell_type <- rownames(df)
df_long <- df %>%
  pivot_longer(cols = -c(cell_type), 
               names_to = "Comparison",
               values_to = 'Num_Degs')
cell_type_order <- rownames(df[order(rowSums(df[,1:3]), decreasing = TRUE),])
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

# Grab colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[nrow(celltype_colors) + 1,] <- c("shared", "gray", "gray")
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
celltype_colors <- celltype_colors[match(cell_type_order, celltype_colors$predicted.celltype.l2),]

# Generate plot
p <- ggplot(df_long, aes(x = Comparison, y = Num_Degs, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = celltype_colors$new_color) +
  theme_Publication_blank() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right",
        text = element_text(size = 24))
p

set_panel_size(p, file = paste0(output_dir, "ADvsHC_inAPOE_uniqueDEGs_celltype_origin.pdf"),
               width = unit(4, "in"), height = unit(7, "in"))

#-------------------------------------------------------------------------------
# 4carriervs33 in Diagnosis # DEGs heatmap
#-------------------------------------------------------------------------------

# Define folders to compare
comparisons <- c("4carriervs33_HC", "4carriervs33_AD")
deg_list <- list()
gene_names_list <- list()

# Define cell types
celltype_list <- c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                   "CD14_Mono", "CD16_Mono",
                   "ASDC", "cDC1", "cDC2", "pDC",
                   "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                   "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                   "NK", "NK_Proliferating", "NK_CD56bright",
                   "ILC", "dnT", "gdT")

for (comparison in comparisons) {
  # Read in intersection
  intersected_df <- read.csv(paste0(input_dir, comparison, "_",
                                    "padj", as.character(padj.thresh),
                                    "_lfc", as.character(lfc.thresh),
                                    "_MAST+edgeR_DEGs.csv"), row.names = 1)
  
  # Create list
  intersected_list <- list()
  for(cell_type in names(intersected_df)) {
    intersected_list[[cell_type]] <- intersected_df[ , cell_type][ (intersected_df[ , cell_type] != "") & (!is.na(intersected_df[ , cell_type])) ]
  }
  
  # Initialize deg vector
  deg_vector <- c()
  gene_names <- c()
  
  for (cell_type in celltype_list) {
    # Add sig genes to list
    if (cell_type %in% names(intersected_list)) {
      deg_vector <- c(deg_vector, length(intersected_list[[cell_type]]))
      gene_names <- union(gene_names, intersected_list[[cell_type]])
    } else {
      deg_vector <- c(deg_vector, NA)
    }
  }

  # Add genes to list
  deg_list[[comparison]] <- deg_vector
  gene_names_list[[comparison]] <- gene_names
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- comparisons
data

# Remove zero rows and all NA rows
data <- data[-which(rowSums(data) == 0),]
data <- data[rowSums(is.na(data)) != ncol(data), ]

# Generate with hclust
pdf(paste0(output_dir, "4carriervs33_numDEGs_heatmap.pdf"),
    width = 5, height = 3)
data_heatmap <- t(data[,1:2])
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         border_color = NA
)
dev.off()

# Make upset for confirmation
pdf(paste0(output_dir, "4carriervs33_numDEGs_upset.pdf"))
print(upset(
  fromList(gene_names_list),
  text.scale = 2,
  sets = names(gene_names_list),
  nsets = length(gene_names_list),
  sets.bar.color = c("gray", "red")
))
dev.off()

#-------------------------------------------------------------------------------
# ADvsHC in APOE # DARs heatmap
#-------------------------------------------------------------------------------

# Define folders to compare
comparisons <- c("ADvsHC_33", "ADvsHC_34", "ADvsHC_44")
deg_list <- list()
gene_names_list <- list()

# Define cell types
celltype_list <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4._T_Cells",
                "CD8._T_Cells", "NK_Cells", "Other_T_Cells", "Other")

for (comparison in comparisons) {
  # Read in intersection
  intersected_df <- read.csv(paste0(da_input_dir, comparison, "_",
                                    "padj", as.character(padj.thresh),
                                    "_lfc", as.character(lfc.thresh),
                                    "_LR+DESeq2_DARs.csv"), row.names = 1)
  
  # Create list
  intersected_list <- list()
  for(cell_type in names(intersected_df)) {
    intersected_list[[cell_type]] <- intersected_df[ , cell_type][ (intersected_df[ , cell_type] != "") & (!is.na(intersected_df[ , cell_type])) ]
  }
  
  # Initialize deg vector
  deg_vector <- c()
  gene_names <- c()
  
  for (cell_type in celltype_list) {
    # Add sig genes to list
    if (cell_type %in% names(intersected_list)) {
      deg_vector <- c(deg_vector, length(intersected_list[[cell_type]]))
      gene_names <- union(gene_names, intersected_list[[cell_type]])
    } else {
      deg_vector <- c(deg_vector, NA)
    }
  }

  # Add genes to list
  deg_list[[comparison]] <- deg_vector
  gene_names_list[[comparison]] <- gene_names
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- c("E3/E3", "E3/E4", "E4/E4")
data

# Remove zero rows and all NA rows
data <- data[-which(rowSums(data) == 0),]
data <- data[rowSums(is.na(data)) != ncol(data), ]

# Generate with hclust
pdf(paste0(output_dir, "ADvsHC_inAPOE_numDARs_heatmap.pdf"),
    width = 5, height = 3)
data_heatmap <- t(data[,1:3])
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         border_color = NA
)
dev.off()

# Make upset for confirmation
pdf(paste0(output_dir, "ADvsHC_inAPOE_numDARs_upset.pdf"))
print(upset(
  fromList(gene_names_list),
  text.scale = 2,
  sets = names(gene_names_list),
  nsets = length(gene_names_list),
  sets.bar.color = c("red", "orange", "yellow")
))
dev.off()

# Identify unique APOE DEGs
unique_gene_names_list <- list()
unique_gene_names_list[["ADvsHC_33"]] <- setdiff(gene_names_list[["ADvsHC_33"]],
                                                 union(gene_names_list[["ADvsHC_34"]], gene_names_list[["ADvsHC_44"]]))
unique_gene_names_list[["ADvsHC_34"]] <- setdiff(gene_names_list[["ADvsHC_34"]],
                                                 union(gene_names_list[["ADvsHC_33"]], gene_names_list[["ADvsHC_44"]]))
unique_gene_names_list[["ADvsHC_44"]] <- setdiff(gene_names_list[["ADvsHC_44"]],
                                                 union(gene_names_list[["ADvsHC_34"]], gene_names_list[["ADvsHC_33"]]))

# Initialize cell type counter
df <- data.frame(row.names = c(celltype_list, "shared"),
                 "ADvsHC_33" = rep(0, length(celltype_list) + 1),
                 "ADvsHC_34" = rep(0, length(celltype_list) + 1),
                 "ADvsHC_44" = rep(0, length(celltype_list) + 1))

for (comparison in comparisons) {
  # Read in intersection
  intersected_df <- read.csv(paste0(da_input_dir, comparison, "_",
                                    "padj", as.character(padj.thresh),
                                    "_lfc", as.character(lfc.thresh),
                                    "_LR+DESeq2_DARs.csv"), row.names = 1)
  
  # Create list
  intersected_list <- list()
  for(cell_type in names(intersected_df)) {
    intersected_list[[cell_type]] <- intersected_df[ , cell_type][ (intersected_df[ , cell_type] != "") & (!is.na(intersected_df[ , cell_type])) ]
  }
  
  # Specify DEGs
  gene_names <- unique_gene_names_list[[comparison]]
  for (gene in gene_names) {
    # initialize result vector
    gene_celltypes <- c()
    for (cell_type in celltype_list[celltype_list %in% names(intersected_list)]) {
      if (gene %in% intersected_list[[cell_type]]) {
        gene_celltypes <- c(gene_celltypes, cell_type)
      }
    }
    if (length(gene_celltypes) == 1) {
      df[gene_celltypes[1], comparison] <- df[gene_celltypes[1], comparison] + 1
    } else {
      df["shared", comparison] <- df["shared", comparison] + 1
    }
  }
}

# Convert to long format
df$cell_type <- rownames(df)
df_long <- df %>%
  pivot_longer(cols = -c(cell_type), 
               names_to = "Comparison",
               values_to = 'Num_Dars')
cell_type_order <- rownames(df[order(rowSums(df[,1:3]), decreasing = TRUE),])
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

# Grab colors
celltype_colors <- read.csv(da_celltype_colors_path)
celltype_colors[nrow(celltype_colors) + 1,] <- c("shared", "gray")
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub("\\+", ".", x)})
celltype_colors <- celltype_colors[match(cell_type_order, celltype_colors$predicted.celltype.l2),]

# Generate plot
p <- ggplot(df_long, aes(x = Comparison, y = Num_Dars, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = celltype_colors$new_color) +
  theme_Publication_blank() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "right",
        text = element_text(size = 24))
p

set_panel_size(p, file = paste0(output_dir, "ADvsHC_inAPOE_uniqueDARs_celltype_origin.pdf"),
               width = unit(4, "in"), height = unit(7, "in"))

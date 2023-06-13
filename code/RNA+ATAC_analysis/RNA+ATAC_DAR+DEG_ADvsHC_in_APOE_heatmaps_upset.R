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
celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/celltype_color_map.csv"
da_celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/broad_celltype_color_map.csv"
input_dir <- "/projects/b1169/projects/AD_APOE/results/de/MAST_edgeR/out_NP_05-19-2023/full_intersection/"
da_input_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/LR_DESeq2/out_NP_05-22-2023/full_intersection/"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_de_comparison/APOE_heatmaps_upset/out_NP_05-23-2023/"
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

# # Get deg to gene ratio
# gene_nums <- sapply(cell_types, grab_gene_num, de_dir = mast_dir, USE.NAMES = TRUE) %>% as.data.frame
# names(gene_nums) <- "num_genes"
# deg_nums <- apply(intersected_df, MARGIN = 2, function(x) {length(x[!is.na(x) & x != ""])}) %>% as.data.frame
# names(deg_nums) <- "num_degs"
# gene_nums <- merge(gene_nums, deg_nums, by = 0, all.y = TRUE)
# gene_nums$ratio <- gene_nums$num_degs / gene_nums$num_genes
# gene_nums <- gene_nums[order(gene_nums$ratio),,drop=FALSE]
# 
# # Grab cell type order
# gene_nums <- gene_nums[which(gene_nums$Row.names %in% names(intersected_list)),]
# celltype_order <- gene_nums$Row.names
# 
# # Read in color map
# color_map <- read.csv(celltype_colors_path)
# color_map$predicted.celltype.l2 <- sapply(color_map$predicted.celltype.l2,
#                                           function(x) {gsub(" ", "_", x)})
# 
# # Change order
# num_deg_order <- gene_nums[order(gene_nums$num_degs), "Row.names"]
# color_map <- color_map[match(celltype_order, color_map$predicted.celltype.l2),]
# 
# # Create upset plot and export
# pdf(file = paste0(output_dir, comparison, "_",
#                   "padj", as.character(padj.thresh),
#                   "_lfc", as.character(lfc.thresh),
#                   "_DEGoverlap_upset.pdf"))
# print(upset(
#   fromList(intersected_list),
#   sets = rev(celltype_order),
#   # text.scale = 1.3,
#   nsets = length(intersected_list),
#   # order.by = "freq",
#   keep.order = TRUE,
#   sets.bar.color = rev(color_map$new_color)
# ))
# dev.off()
# 
# # Create heatmap of gene numbers and ratios
# rownames(gene_nums) <- gene_nums$Row.names
# gene_nums <- gene_nums[,"ratio",drop=FALSE]
# 
# pdf(file = paste0(output_base_dir, comparison, "_",
#                   "padj", as.character(padj.thresh),
#                   "_lfc", as.character(lfc.thresh),
#                   "_DEGratio_heatmap.pdf"))
# pheatmap(gene_nums,
#          cluster_rows = FALSE, 
#          cluster_cols = TRUE,
#          display_numbers = TRUE,
#          number_format = "%.4f",
#          cellwidth = 40,
#          color = colorRampPalette(c("white", "tomato2"))(100)
# )
# dev.off()

# #-------------------------------------------------------------------------------
# # Get number of degs for each celltype (inter APOE)
# 
# # Define folders to compare
# comparisons <- c("apoe34vs33", "apoe44vs33", "apoe44vs34")
# deg_list <- list()
# 
# # Get celltype list
# celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
#                         function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- paste0(de_base_dir, comparison, "_ad", "/out_NP_11-10-2022_covarSex/")
# 
#   # Make list of csvs
#   csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
# 
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   # Specify APOE label
#   if (comparison == "apoe34vs33") {
#     label <- "34_vs_33"
#   } else if (comparison == "apoe44vs33") {
#     label <- "44_vs_33"
#   } else { label <- "44_vs_34"}
# 
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_ADonly_degs.csv"))
# 
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
# 
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, length(sig_genes$X))
#     }, error = function(e) {
#       NULL
#     })
#   }
# 
#   # Add genes to list
#   # deg_list[[comparison]] <- unique(deg_vector)
#   deg_list[[comparison]] <- deg_vector
# }
# 
# # Convert to dataframe
# data <- as.data.frame(deg_list)
# rownames(data) <- celltype_list
# colnames(data) <- comparisons
# ad_data <- data
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- paste0(de_base_dir, comparison, "_hc", "/out_NP_11-10-2022_covarSex/")
#   
#   # Make list of csvs
#   csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   # Specify APOE label
#   if (comparison == "apoe34vs33") {
#     label <- "34_vs_33"
#   } else if (comparison == "apoe44vs33") {
#     label <- "44_vs_33"
#   } else { label <- "44_vs_34"}
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_HConly_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, length(sig_genes$X))
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   # deg_list[[comparison]] <- unique(deg_vector)
#   deg_list[[comparison]] <- deg_vector
# }
# 
# # Convert to dataframe
# data <- as.data.frame(deg_list)
# rownames(data) <- celltype_list
# colnames(data) <- comparisons
# hc_data <- data
# 
# # Find difference
# diff <- ad_data / hc_data
# 
# # Calculate LFC
# diff <- log2(diff)
# 
# # Remove nonfinite rows
# diff <- diff[is.finite(rowSums(diff)),]
# 
# # Add celltype column
# diff$celltype <- rownames(diff)
# 
# # Convert to long format
# diff_long <- pivot_longer(
#   diff,
#   cols = `apoe34vs33`:`apoe44vs34`,
#   names_to = "APOE",
#   values_to = "Log2Fold_Change_ADoverHC"
# )
# 
# # Generate plot
# p <- ggplot(diff_long, aes(x = celltype, y = APOE,
#                            fill = Log2Fold_Change_ADoverHC)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   theme_Publication_blank() +
#   scale_x_discrete(position = 'top', expand = c(0,0)) +
#   scale_y_discrete(limits = rev, expand = c(0,0)) +
#   theme(legend.position = "right",
#         legend.direction="vertical",
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
# p
# 
# # Export plot
# set_panel_size(p, file = paste0(output_dir, "pbmcsRNA_APOEcomp_FC-OfDegsADvsHC.pdf"),
#                width = unit(8, "in"), height = unit(2, "in"))
# 
# # Generate with hclust
# data_heatmap <- t(diff[,1:3])
# colfunc <- colorRampPalette(c("blue", "white", "red"))
# pheatmap(as.matrix(data_heatmap),
#          cluster_rows = FALSE,
#          color = colfunc(100),
#          breaks=seq(-5, 5, length.out=101),
#          border_color = NA,
#          filename = paste0(output_dir, "pbmcsRNA_APOEcomp_FC-OfDegsADvsHC_clust.pdf"),
#          width = 5,
#          height = 3
# )
# 
# #-------------------------------------------------------------------------------
# # Get number of degs for each celltype (inter APOE, just num)
# 
# # Define folders to compare
# comparisons <- c("apoe34vs33", "apoe44vs33", "apoe44vs34")
# deg_list <- list()
# 
# # Get celltype list
# celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
#                         function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- paste0(de_base_dir, comparison, "_ad", "/out_NP_11-10-2022_covarSex/")
#   
#   # Make list of csvs
#   csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   # Specify APOE label
#   if (comparison == "apoe34vs33") {
#     label <- "34_vs_33"
#   } else if (comparison == "apoe44vs33") {
#     label <- "44_vs_33"
#   } else { label <- "44_vs_34"}
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_ADonly_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, length(sig_genes$X))
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   # deg_list[[comparison]] <- unique(deg_vector)
#   deg_list[[comparison]] <- deg_vector
# }
# ad_deg_list <- deg_list
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- paste0(de_base_dir, comparison, "_hc", "/out_NP_11-10-2022_covarSex/")
#   
#   # Make list of csvs
#   csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   # Specify APOE label
#   if (comparison == "apoe34vs33") {
#     label <- "34_vs_33"
#   } else if (comparison == "apoe44vs33") {
#     label <- "44_vs_33"
#   } else { label <- "44_vs_34"}
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_HConly_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, length(sig_genes$X))
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   # deg_list[[comparison]] <- unique(deg_vector)
#   deg_list[[comparison]] <- deg_vector
# }
# hc_deg_list <- deg_list
# 
# # Convert to dataframe
# data_ad <- as.data.frame(ad_deg_list)
# data_hc <- as.data.frame(hc_deg_list)
# colnames(data_hc) <- paste0("HC_", comparisons)
# colnames(data_ad) <- paste0("AD_", comparisons)
# data <- cbind(data_hc, data_ad)
# data <- data[,c(1,4,3,6,2,5)]
# rownames(data) <- celltype_list
# 
# # Generate with hclust
# data_heatmap <- t(data)
# colfunc <- colorRampPalette(c("white", "red"))
# pheatmap(as.matrix(data_heatmap),
#          cluster_rows = FALSE,
#          color = colfunc(100),
#          # breaks=seq(-5, 5, length.out=101),
#          border_color = NA,
#          filename = paste0(output_dir, "pbmcsRNA_APOEcomp_AD_HC_numDegs_clust.pdf"),
#          width = 5,
#          height = 3
# )
# 
# #-------------------------------------------------------------------------------
# # Get number of degs for each celltype (inter APOE, just num, clonal only)
# 
# # Define folders to compare
# comparisons <- c("apoe34_vs_33", "apoe44_vs_33", "apoe44_vs_34")
# deg_list <- list()
# 
# # Get celltype list
# celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
#                         function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- "/projects/b1169/projects/AD_APOE/results/tcr/clonality/out_2022_12_14_AR/"
# 
#   # Initialize deg vector
#   deg_vector <- c()
# 
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, "Alzheimers_Disease_", cell_type, "_clonal_", comparison, "_degs.csv"))
# 
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
# 
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, nrow(sig_genes))
#     }, error = function(e) {
#       NULL
#     })
#   }
# 
#   # Add genes to list
#   deg_list[[comparison]] <- deg_vector
# }
# ad_deg_list <- deg_list
# 
# cell_types <- c()
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- "/projects/b1169/projects/AD_APOE/results/tcr/clonality/out_2022_12_14_AR/"
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, "Healthy_Control_", cell_type, "_clonal_", comparison, "_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, nrow(sig_genes))
#       
#       cell_types <- c(cell_types, cell_type)
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   deg_list[[comparison]] <- deg_vector
# }
# hc_deg_list <- deg_list
# 
# # Convert to dataframe
# data_ad <- as.data.frame(ad_deg_list)
# data_hc <- as.data.frame(hc_deg_list)
# colnames(data_hc) <- paste0("HC_", comparisons)
# colnames(data_ad) <- paste0("AD_", comparisons)
# data <- cbind(data_hc, data_ad)
# data <- data[,c(1,4,3,6,2,5)]
# rownames(data) <- cell_types[1:9]
# data <- data[rowSums(data) > 0,]
#  
# # Generate with hclust
# data_heatmap <- t(data)
# colfunc <- colorRampPalette(c("white", "red"))
# pheatmap(as.matrix(data_heatmap),
#          cluster_rows = TRUE,
#          color = colfunc(100),
#          # breaks=seq(-5, 5, length.out=101),
#          border_color = NA,
#          filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_numDegs_clust_RowClust.pdf"),
#          width = 5,
#          height = 3
# )
# 
# #-------------------------------------------------------------------------------
# # Get number of degs for each celltype (inter APOE, lfc, clonal only)
# 
# # Define folders to compare
# comparisons <- c("apoe34_vs_33", "apoe44_vs_33", "apoe44_vs_34")
# deg_list <- list()
# 
# # Get celltype list
# celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
#                         function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()
# 
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- "/projects/b1169/projects/AD_APOE/results/tcr/clonality/out_2022_12_14_AR/"
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, "Alzheimers_Disease_", cell_type, "_clonal_", comparison, "_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, nrow(sig_genes))
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   deg_list[[comparison]] <- deg_vector
# }
# ad_deg_list <- deg_list
# 
# cell_types <- c()
# for (comparison in comparisons) {
#   # Define de dir
#   # comparison <- comparisons[1]
#   de_dir <- "/projects/b1169/projects/AD_APOE/results/tcr/clonality/out_2022_12_14_AR/"
#   
#   # Initialize deg vector
#   deg_vector <- c()
#   
#   for (cell_type in celltype_list) {
#     # Load in degs
#     tryCatch({
#       # Load in degs
#       degs <- read.csv(paste0(de_dir, "Healthy_Control_", cell_type, "_clonal_", comparison, "_degs.csv"))
#       
#       # Identify sig genes
#       sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
#       
#       # Add sig genes to list
#       deg_vector <- c(deg_vector, nrow(sig_genes))
#       
#       cell_types <- c(cell_types, cell_type)
#     }, error = function(e) {
#       NULL
#     })
#   }
#   
#   # Add genes to list
#   deg_list[[comparison]] <- deg_vector
# }
# hc_deg_list <- deg_list
# 
# # # Convert to dataframe
# # data_ad <- as.data.frame(ad_deg_list)
# # data_hc <- as.data.frame(hc_deg_list)
# # colnames(data_hc) <- paste0("HC_", comparisons)
# # colnames(data_ad) <- paste0("AD_", comparisons)
# # data <- cbind(data_hc, data_ad)
# # data <- data[,c(1,4,3,6,2,5)]
# # rownames(data) <- cell_types[1:9]
# # data <- data[rowSums(data) > 0,]
# 
# # Convert to dataframe
# data_ad <- as.data.frame(ad_deg_list)
# data_hc <- as.data.frame(hc_deg_list)
# colnames(data_hc) <- paste0("HC_", comparisons)
# colnames(data_ad) <- paste0("AD_", comparisons)
# rownames(data_hc) <- cell_types[1:9]
# rownames(data_ad) <- cell_types[1:9]
# 
# # Find difference
# diff <- data_ad / data_hc
# 
# # Calculate LFC
# diff <- log2(diff)
# 
# # Remove nonfinite rows
# diff <- diff[is.finite(rowSums(diff)),]
# 
# # Generate with hclust
# data_heatmap <- t(data)
# colfunc <- colorRampPalette(c("white", "red"))
# # colfunc <- colorRampPalette(c("white", "red", "blue"))
# pheatmap(as.matrix(data_heatmap),
#          cluster_rows = FALSE,
#          color = colfunc(100),
#          # breaks=seq(-5, 5, length.out=101),
#          border_color = NA,
#          filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_numDegs_clust_HCtop.pdf"),
#          # filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_LFC_clust_RowClust.pdf"),
#          width = 5,
#          height = 3
# )
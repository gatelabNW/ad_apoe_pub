# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-18-2022
# Written by: Abhi Ramakrishnan
# Summary: Make umaps colored by T and B cell clonality of scRNA seq data
# 
# ------------------------------------------------------------------------------

# Load packages
suppressMessages({
  library("plyr")
  library("dplyr")
  library("tidyverse")
  library("Seurat")
  library("ggplot2")
  library("ggrastr")
})

# Create output directory
output_dir <- "/path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# source helper functions
source("path/to/helper_functions.R")

# ------------------------------------------------------------------------------
# UMAP annotated by T-cell clonality
# ------------------------------------------------------------------------------
# Load seurat object with dimension reduction "umap"
s <- "/path/to/object/"
load(s)

# Get UMAP coordinates
umap_coords <- data.frame(Embeddings(s, reduction = "umap"))

# Load seurat object with T cell data
s <- "/path/to/object/"
load(s)

# Get clonality information
clonal_df <- s[["clonal_T"]]

# Check that barcodes are in same order 
all(rownames(umap_coords) == rownames(clonal_df))

# Bind umap coords and metadata
clonal_df <- cbind(umap_coords, clonal_df)

# Eliminate cells that do not have clonality information
clonal_df <- clonal_df %>% dplyr::filter(clonal_T != "NA")

# Plot umap

p <- ggplot(data = clonal_df, 
            aes(x=UMAP_1, y=UMAP_2, color=clonal_T)) +
  geom_point(size=.2) +
  scale_color_manual(values= c("#6a51a3", "#fb6a4a")) +
  theme_Publication_blank() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p
set_panel_size(p, file= paste0(output_dir, "umap_T_clonal_v_nonclonal.pdf"))


# ------------------------------------------------------------------------------
# UMAP annotated by B-cell clonality
# ------------------------------------------------------------------------------
# Load seurat object with dimension reduction "umap"
s <- "/path/to/object/"
# Get UMAP coordinates
umap_coords <- data.frame(Embeddings(s, reduction = "umap"))

# Load seurat object with B cell data
s <- "/path/to/object/"
load(s)

# Get clonality information and cell subtypes
clonal_df <- s[[c("clonal_B", "predicted.celltype.l2")]]

# Check that barcodes are in same order 
all(rownames(umap_coords) == rownames(clonal_df))
# Bind umap coords and metadata
clonal_df <- cbind(umap_coords, clonal_df)

# Eliminate cells that do not have clonality information
clonal_df <- clonal_df %>% dplyr::filter(clonal_B != "NA")

# 1. Plot umap colored by clonality

p <- ggplot(data = clonal_df, 
            aes(x=UMAP_1, y=UMAP_2, color=clonal_B)) +
  geom_point(size=.2) +
  scale_color_manual(values= c("#6a51a3", "#fb6a4a")) +
  theme_Publication_blank() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p
set_panel_size(p, file= paste0(output_dir, "umap_B_clonal_v_nonclonal.pdf"))

# 2. Plot umap colored by B cell subtype

p <- ggplot(data = clonal_df, 
            aes(x=UMAP_1, y=UMAP_2, color=predicted.celltype.l2)) +
  geom_point(size=.2) +
  scale_color_manual(values= c("#f58383", "#fbcdcd", "#Eb0707", "#a50505")) +
  theme_Publication_blank() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p
set_panel_size(p, file= paste0(output_dir, "umap_B_cell_subtype.pdf"))





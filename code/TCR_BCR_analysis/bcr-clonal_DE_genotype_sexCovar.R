# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-10-2023
# Written by: Abhi Ramakrishnan, modified from Natalie Piehl's DE scripts
# Summary: Perform differential expression analysis of clonal B-cells by diagnosis and genotype for AD-APOE project
# Use Sex as covariate
# ------------------------------------------------------------------------------
# Create output directory
output_dir <- "path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# Source Natalie's helper functions
source("path/to/helper_functions.R")

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("doMC")
  library("UpSetR")
})

# Load seurat object
s <- "path/to/object/"
load(s)

# ------------------------------------------------------------------------------
# 1. Subset object by clonality and APOE genotype
# ------------------------------------------------------------------------------
# List all B cell clusters which have enough cells for DE analysis
b_cell_clusters <- c(
  "B naive",
  "B intermediate",
  "B memory",
  "Plasmablast"
)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)
# Set identity used for DE
s <- SetIdent(s, value = "Diagnosis")
# Subset Seurat object by clonality
object_clonal <- subset(s, clonal_B== "C")
# Further subset object by APOE genotype
object_clonal_33 <- subset(object_clonal, APOE_genotype == "E3/E3")
object_clonal_34 <- subset(object_clonal, APOE_genotype == "E3/E4")
object_clonal_44 <- subset(object_clonal, APOE_genotype == "E4/E4")
# # Subset clonal object by cell type and verify subsetting works
print(table(object_clonal_33$clonal_B, object_clonal_33$predicted.celltype.l2))

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25
# ------------------------------------------------------------------------------
# 2. Perform DE analysis of AD vs HC for all clonal B cells by genotype
# ------------------------------------------------------------------------------
# APOE genotype E3/E3
x <- FindMarkers(object = object_clonal_33,
                 ident.1 = "Alzheimers Disease",
                 ident.2= "Healthy Control",
                 logfc.threshold = -Inf,
                 test.use = "MAST",
                 min.pct = 0.1,
                 assay= 'RNA',
                 latent.vars = "Sex"
)
# Remove ribosomal, mitochondrial (leaving HLA genes present for now)
if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
  x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
}
# Run Benjamini-Hochberg adjustment
x$BH <- p.adjust(x$p_val, method = "BH")
print(head(x))

# Write out results
write.csv(x, paste0(output_dir, "all_AD_v_HC_E3-E3_clonal_degs.csv"))
# APOE genotype E3/E4
x <- FindMarkers(object = object_clonal_34,
                 ident.1 = "Alzheimers Disease",
                 ident.2= "Healthy Control",
                 logfc.threshold = -Inf,
                 test.use = "MAST",
                 min.pct = 0.1,
                 assay= 'RNA',
                 latent.vars = "Sex"
)
# Remove ribosomal, mitochondrial (leaving HLA genes present for now)
if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
  x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
}
# Run Benjamini-Hochberg adjustment
x$BH <- p.adjust(x$p_val, method = "BH")
print(head(x))

# Write out results
write.csv(x, paste0(output_dir, "all_AD_v_HC_E3-E4_clonal_degs.csv"))
# APOE genotype E4/E4
x <- FindMarkers(object = object_clonal_44,
                 ident.1 = "Alzheimers Disease",
                 ident.2= "Healthy Control",
                 logfc.threshold = -Inf,
                 test.use = "MAST",
                 min.pct = 0.1,
                 assay= 'RNA',
                 latent.vars = "Sex"
)
# Remove ribosomal, mitochondrial (leaving HLA genes present for now)
if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
  x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
}
# Run Benjamini-Hochberg adjustment
x$BH <- p.adjust(x$p_val, method = "BH")
print(head(x))

# Write out results
write.csv(x, paste0(output_dir, "all_AD_v_HC_E4-E4_clonal_degs.csv"))
# ------------------------------------------------------------------------------
# 3. Run DE for each cell type for APOE genotype E3/E3
# ------------------------------------------------------------------------------
# Redefine B cell clusters which have enough cells for DE analysis
b_cell_clusters <- c(
  "B naive",
  "B intermediate",
  "B memory"
)
for (cluster in b_cell_clusters) {
  x <- FindMarkers(object = subset(object_clonal_33, predicted.celltype.l2 == cluster),
                   ident.1 = "Alzheimers Disease",
                   ident.2= "Healthy Control",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = "Sex"
  )
  assign(paste0("de_AD_v_HC_E3-E3_clonal_", cluster), x)
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))
  
  # Write out results
  cluster_name <- gsub(" ", "_", cluster)
  write.csv(x, paste0(output_dir, cluster_name, "_AD_v_HC_E3-E3_clonal_degs.csv"))
  
  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E3/E3 Clonal AD vs HC in ", cluster),
               file = paste0(output_dir, cluster, "_AD_v_HC_E3-E3_clonal_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  
}
# ------------------------------------------------------------------------------
# 4. Run DE for each cell type for APOE genotype E3/E4
# ------------------------------------------------------------------------------
# Redefine B cell clusters which have enough cells for DE analysis
b_cell_clusters <- c(
  "B naive",
  "B intermediate",
  "B memory"
)
for (cluster in b_cell_clusters) {
  x <- FindMarkers(object = subset(object_clonal_34, predicted.celltype.l2 == cluster),
                   ident.1 = "Alzheimers Disease",
                   ident.2= "Healthy Control",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = "Sex"
  )
  assign(paste0("de_AD_v_HC_E3-E4_clonal_", cluster), x)
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))
  
  # Write out results
  cluster_name <- gsub(" ", "_", cluster)
  write.csv(x, paste0(output_dir, cluster_name, "_AD_v_HC_E3-E4_clonal_degs.csv"))
  
  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E3/E4 Clonal AD vs HC in ", cluster),
               file = paste0(output_dir, cluster, "_AD_v_HC_E3-E4_clonal_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  
}
# ------------------------------------------------------------------------------
# 5. Run DE for each cell type for APOE genotype E4/E4
# ------------------------------------------------------------------------------
# Redefine B cell clusters which have enough cells for DE analysis
b_cell_clusters <- c(
  "B naive",
  "B intermediate",
  "B memory"
)
for (cluster in b_cell_clusters) {
  x <- FindMarkers(object = subset(object_clonal_44, predicted.celltype.l2 == cluster),
                   ident.1 = "Alzheimers Disease",
                   ident.2= "Healthy Control",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = "Sex"
  )
  assign(paste0("de_AD_v_HC_E4-E4_clonal_", cluster), x)
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))
  
  # Write out results
  cluster_name <- gsub(" ", "_", cluster)
  write.csv(x, paste0(output_dir, cluster_name, "_AD_v_HC_E4-E4_clonal_degs.csv"))
  
  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E4/E4 Clonal AD vs HC in ", cluster),
               file = paste0(output_dir, cluster, "_AD_v_HC_E4-E4_clonal_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  
}
# ------------------------------------------------------------------------------
# 6. Create upset plots for all DEGs by genotype (modified from Natalie Piehl's upset plot script)
# ------------------------------------------------------------------------------
sig_genes_ls <- list()
apoe_list <- list("E3-E3",
                  "E3-E4",
                  "E4-E4")
# Create list with sig genes for each cell type
for (genotype in apoe_list) {
  print(genotype)
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, "all_AD_v_HC_", genotype, "_clonal_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[genotype]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get colors
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_genotype.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = c("orange", "yellow", "red")
)
dev.off()
# ------------------------------------------------------------------------------
# 7. Create upset plots for cell type specific DEGs by genotype (modified from Natalie Piehl's upset plot script)
# ------------------------------------------------------------------------------
# APOE E3/E3
sig_genes_ls <- list()
# Create list with sig genes for each cell type
for (cluster in b_cell_clusters) {
  cluster_name <- gsub(" ", "_", cluster)
  print(cluster_name)
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, cluster_name, "_AD_v_HC_E3-E3", "_clonal_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[cluster]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get colors
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype_E3-E3.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()

# APOE E3/E4
sig_genes_ls <- list()
# Create list with sig genes for each cell type
for (cluster in b_cell_clusters) {
  cluster_name <- gsub(" ", "_", cluster)
  print(cluster_name)
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, cluster_name, "_AD_v_HC_E3-E4", "_clonal_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[cluster]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get colors
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype_E3-E4.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()

# APOE E4/E4
sig_genes_ls <- list()
# Create list with sig genes for each cell type
for (cluster in b_cell_clusters) {
  cluster_name <- gsub(" ", "_", cluster)
  print(cluster_name)
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, cluster_name, "_AD_v_HC_E4-E4", "_clonal_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[cluster]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get colors
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype_E4-E4.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()
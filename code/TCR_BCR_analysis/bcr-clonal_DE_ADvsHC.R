# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-09-2022
# Written by: Abhi Ramakrishnan, modified from Natalie Piehl's DE scripts
# Summary: Perform differential expression analysis of clonal B-cells by diagnosis for AD-APOE project
# 
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
# 1. Subset seurat object by clonal B cell groups and run differential expression by diagnosis
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
# verify subsetting works
print(unique(object_clonal@meta.data$predicted.celltype.l2))

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Run DE for clonal HC vs clonal AD cells for each B cell subtype
for (cluster in b_cell_clusters) {
  x <- FindMarkers(object = subset(object_clonal, predicted.celltype.l2 == cluster),
                   ident.1 = "Alzheimers Disease",
                   ident.2= "Healthy Control",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = c("Sex", "APOE_genotype")
                   )
  assign(paste0("de_diagnosis_clonal_", cluster), x)
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))

  # Write out results
  write.csv(x, paste0(output_dir, cluster, "_clonal_AD_vs_HC_degs.csv"))

  # Create volcano plot
  volcano_plot(x, title = paste0("Clonal AD vs HC in ", cluster),
               file = paste0(output_dir, cluster, "_clonal_AD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(cluster, " done"))

}

# ------------------------------------------------------------------------------
# 2. Subset seurat object by nonclonal Tcell groups and run differential expression by diagnosis
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
object_nonclonal <- subset(s, clonal_B== "NC")
print(unique(object_nonclonal@meta.data$predicted.celltype.l2))

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Run DE for clonal HC vs clonal AD cells for each B cell subtype
for (cluster in b_cell_clusters) {
  x <- FindMarkers(object = subset(object_nonclonal, predicted.celltype.l2 == cluster),
                   ident.1 = "Alzheimers Disease",
                   ident.2= "Healthy Control",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = c("Sex", "APOE_genotype")
  )
  assign(paste0("de_diagnosis_nonclonal_", cluster), x)
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))
  
  # Write out results
  write.csv(x, paste0(output_dir, cluster, "_nonclonal_AD_vs_HC_degs.csv"))
  
  # Create volcano plot
  volcano_plot(x, title = paste0("Nonclonal AD vs HC in ", cluster),
               file = paste0(output_dir, cluster, "_nonclonal_AD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(cluster, " done"))
  
}

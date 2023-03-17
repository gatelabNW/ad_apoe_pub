# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-26-2023
# Written by: Abhi Ramakrishnan, modified from Natalie Piehl's DE scripts
# Summary: Perform differential expression analysis of clonal B-cells by diagnosis for AD-APOE project
# Sex covariate
# ------------------------------------------------------------------------------
# Create output directory
output_dir <- "path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# Source Natalie's helper functions
source("/path/to/helper_functions.R")

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("UpSetR")
})
# Load seurat object
s <- "path/to/object/"
load(s)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set identity used for DE
s <- SetIdent(s, value = "APOE_genotype")

# Subset Seurat object by B-cell clonality
object_clonal <- subset(s, clonal_B== "C")

# Verify subsetting
print(table(object_clonal$clonal_B, object_clonal$predicted.celltype.l2))

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# List diagnoses
diag_list <- c("Healthy Control", "Alzheimers Disease")

# ------------------------------------------------------------------------------
# 1. Run DE for each genotype comparison within each diagnosis subset
# ------------------------------------------------------------------------------
# Run DE for APOE 3/4 vs 3/3
for (diag in diag_list) {
  print(diag)
  x <- FindMarkers(object = subset(object_clonal, Diagnosis==diag),
                   ident.1 = "E3/E4",
                   ident.2= "E3/E3",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = 'Sex'
  )
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))

  # Write out results
  diag <- gsub(" ", "_", diag)
  write.csv(x, paste0(output_dir, diag, "_clonal_apoe34_vs_33_degs.csv"))

  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E3/E4 vs E3/E3 in ", diag),
               file = paste0(output_dir, diag, "_clonal_apoe34_vs_33_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(diag, " apoe 3/4 vs 3/3 is done!"))
}

# Run DE for APOE 4/4 vs 3/4

for (diag in diag_list) {
  print(diag)
  x <- FindMarkers(object = subset(object_clonal, Diagnosis==diag),
                   ident.1 = "E4/E4",
                   ident.2= "E3/E4",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = 'Sex'
  )
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))

  # Write out results
  diag <- gsub(" ", "_", diag)
  write.csv(x, paste0(output_dir, diag, "_clonal_apoe44_vs_34_degs.csv"))

  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E4/E4 vs E3/E4 in ", diag),
               file = paste0(output_dir, diag, "_clonal_apoe44_vs_34_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(diag, " apoe 4/4 vs 3/4 is done!"))
}

# Run DE for APOE 4/4 vs 3/3

for (diag in diag_list) {
  print(diag)
  x <- FindMarkers(object = subset(object_clonal, Diagnosis==diag),
                   ident.1 = "E4/E4",
                   ident.2= "E3/E3",
                   logfc.threshold = -Inf,
                   test.use = "MAST",
                   min.pct = 0.1,
                   assay= 'RNA',
                   latent.vars = 'Sex'
  )
  # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
  if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
    x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
  }
  # Run Benjamini-Hochberg adjustment
  x$BH <- p.adjust(x$p_val, method = "BH")
  print(head(x))

  # Write out results
  diag <- gsub(" ", "_", diag)
  write.csv(x, paste0(output_dir, diag, "_clonal_apoe44_vs_33_degs.csv"))

  # Create volcano plot
  volcano_plot(x, title = paste0("APOE E4/E4 vs E3/E3 in ", diag),
               file = paste0(output_dir, diag, "_clonal_apoe44_vs_33_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(diag, " apoe 4/4 vs 3/3 is done!"))
}

# ------------------------------------------------------------------------------
# 2. Make upset plots for each diagnosis group to quantify num DEGs per genotype comparison
# ------------------------------------------------------------------------------
# List comparisons
comparisons <- c("apoe34_vs_33",
                 "apoe44_vs_34",
                 "apoe44_vs_33")
colors <- c("#FF90D9",
            "#FFD065",
            "#E769FF")
apoe_comparisons <- data.frame(cbind(comparisons, colors))
# UPset plot for HC
# Create list with sig genes for each cell type
sig_genes_ls <- list()
for (comparison in apoe_comparisons$comparisons) {
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir,"Healthy_Control_clonal_",comparison, "_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(comparison)
    print(length(sig_genes[,1]))
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[comparison]] <- sig_genes$X
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

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

apoe_comparisons <- apoe_comparisons[match(num_degs[,1], apoe_comparisons$comparisons),]
apoe_comparisons <- apoe_comparisons[which(apoe_comparisons$comparisons %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "upset_HC_sex_covar_genotype_comparisons_p0.01_lfc0.125.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = apoe_comparisons$colors
)
dev.off()

# UPset plot for AD
# List comparisons
comparisons <- c("apoe34_vs_33",
                 "apoe44_vs_34",
                 "apoe44_vs_33")
colors <- c("#FF90D9",
            "#FFD065",
            "#E769FF")
apoe_comparisons <- data.frame(cbind(comparisons, colors))
# Create list with sig genes for each cell type
sig_genes_ls <- list()
for (comparison in apoe_comparisons$comparisons) {
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir,"Alzheimers_Disease_clonal_",comparison, "_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(comparison)
    print(length(sig_genes[,1]))
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[comparison]] <- sig_genes$X
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

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

apoe_comparisons <- apoe_comparisons[match(num_degs[,1], apoe_comparisons$comparisons),]
apoe_comparisons <- apoe_comparisons[which(apoe_comparisons$comparisons %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "upset_AD_sex_covar_genotype_comparisons_p0.01_lfc0.125.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = apoe_comparisons$colors
)
dev.off()

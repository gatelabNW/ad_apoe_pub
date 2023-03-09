# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-30-2022
# Written by: Natalie Piehl
# Summary: Make upsets from DARs
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("UpSetR")
})

# Organize inputs
da_base_dir <- "/path/to/da/results"
celltype_colors_path <- "/path/to/celltype_color_map.csv"
output_dir <- "/path/to/output/folder"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Define comparison to look at
comparisons <- c("ADvsHC",
                 "ADvsHC_33", "ADvsHC_34", "ADvsHC_44",
                 "34vs33_HC", "34vs33_AD",
                 "44vs33_HC", "44vs33_AD",
                 "44vs34_HC", "44vs34_AD")
comparison <- comparisons[10]

#------------------------------------------------------------------------------
# Create Upset plot
#-------------------------------------------------------------------------------

# Initialize sig gene list
sig_genes_ls <- list()

# Get cell types
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2, function(x) {gsub(" ", "_", x)})
cell_types <- celltype_colors$predicted.celltype.l2

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  print(cell_type)
  # Load in degs
  tryCatch({
    # Read in peaks
    da_peaks <- read.csv(paste0(da_base_dir, comparison, "/", cell_type, "_", comparison, "_dars.csv"))
    
    # Identify sig regions
    da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & abs(da_peaks$avg_log2FC) > lfc.thresh),]
    
    # Extract DARs
    sig_genes <- da_peaks$sitename

    # Add sig genes to list
    sig_genes_ls[[cell_type]] <- sig_genes
  }, error = function(e) {
    message(e)
  })
}

# Reserve cell types with at least 1 DAR
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]
celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, comparison, "_upset_DARs.pdf"))
print(upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  show.numbers = FALSE,
  sets.bar.color = celltype_colors$new_color,
  text.scale = 2
))
dev.off()
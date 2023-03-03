# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-20-2022
# Written by: Natalie Piehl
# Summary: Make upset from DARs
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("UpSetR")
})

# DARs or genes
dars <- FALSE

# Organize inputs
dar_dir <- "/projects/b1169/projects/AD_APOE/results/da/volcano/out_NP_01-10-2023_full_lfc125/"
celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/celltype_color_map.csv"
output_dir <- "/projects/b1169/projects/AD_APOE/results/da/advshc_apoe_upset/out_NP_01-10-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.125

#------------------------------------------------------------------------------
# Create Upset plot (ADvsHC in APOE)

# Get colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2,
                                                function(x) {gsub(" ", "_", x)})
cell_types <- celltype_colors$predicted.celltype.l2

# Comparisons
comparisons <- c("ADvsHC_33", "ADvsHC_34", "ADvsHC_44")

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  # Initialize sig gene list
  sig_genes_ls <- list()
  tryCatch({
    for (comparison in comparisons) {
      # Load in degs
      tryCatch({
        # Read in peaks
        da_peaks <- read.csv(paste0(dar_dir, comparison, "/", cell_type, "_", comparison, "_annotated_peaks.csv"))
        
        # Identify sig regions
        da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & abs(da_peaks$avg_log2FC) > lfc.thresh),]
        
        if (dars) {
          # Extract regions
          sig_genes <- da_peaks$sitename
        } else {
          # Extract gene names
          sig_genes <- unique(da_peaks$nearestGene)[!is.na(unique(da_peaks$nearestGene))]
        }
        
        # Add sig genes to list
        sig_genes_ls[[comparison]] <- sig_genes
      }, error = function(e) {
        message(e)
      })
    }
    # sig_genes_ls[[cell_type]] <- cell_type_ls
    # Reserve comparisons with at least 1 DAR
    sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]
    
    if (length(sig_genes_ls) == 0) {
      next
    }
    
    # Find number of genes in each set
    num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
    for (key in names(sig_genes_ls)) {
      num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
    }
    
    # Define color map
    color_map <- data.frame(analysis = comparisons,
                            color = c("gold", "orange", "red"))
    
    # Get order of sets for coloring
    num_degs <- num_degs[order(-num_degs[, 2]),]
    color_map <- color_map[match(num_degs[,1], color_map$analysis),]
    color_map <- color_map[ which(color_map$analysis %in% names(sig_genes_ls)),]
    
    # Define dars tag
    if (dars) {
      tag <- "DARs"
    } else {
      tag <- "DAR_nearest_genes"
    }
    
    pdf(file = paste0(output_dir, cell_type, "_ADvsHC_APOE_upset_", tag, ".pdf"))
    print(upset(
      fromList(sig_genes_ls),
      nsets = length(sig_genes_ls),
      order.by = "freq",
      sets.bar.color = color_map$color,
      text.scale = 1.5
    ))
    dev.off()
  }, error = function(e) {
    message(paste(cell_type, e))
  }
  )
}

#------------------------------------------------------------------------------
# Create Upset plot (APOEvsAPOE in AD)

# Comparisons
comparisons <- c("34vs33_AD", "44vs33_AD", "44vs34_AD")

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  # Initialize sig gene list
  sig_genes_ls <- list()
  tryCatch({
    for (comparison in comparisons) {
      # Load in degs
      tryCatch({
        # Read in peaks
        da_peaks <- read.csv(paste0(dar_dir, comparison, "/", cell_type, "_", comparison, "_annotated_peaks.csv"))
        
        # Identify sig regions
        da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & abs(da_peaks$avg_log2FC) > lfc.thresh),]
        
        if (dars) {
          # Extract regions
          sig_genes <- da_peaks$sitename
        } else {
          # Extract gene names
          sig_genes <- unique(da_peaks$nearestGene)[!is.na(unique(da_peaks$nearestGene))]
        }
        
        # Add sig genes to list
        sig_genes_ls[[comparison]] <- sig_genes
      }, error = function(e) {
        message(e)
      })
    }
    # sig_genes_ls[[cell_type]] <- cell_type_ls
    # Reserve comparisons with at least 1 DAR
    sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]
    
    if (length(sig_genes_ls) == 0) {
      next
    }
    
    # Find number of genes in each set
    num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
    for (key in names(sig_genes_ls)) {
      num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
    }
    
    # Define color map
    color_map <- data.frame(analysis = comparisons,
                            color = c("aquamarine", "cyan", "dodgerblue"))
    
    # Get order of sets for coloring
    num_degs <- num_degs[order(-num_degs[, 2]),]
    color_map <- color_map[match(num_degs[,1], color_map$analysis),]
    color_map <- color_map[ which(color_map$analysis %in% names(sig_genes_ls)),]
    
    # Define dars tag
    if (dars) {
      tag <- "DARs"
    } else {
      tag <- "DAR_nearest_genes"
    }
    
    pdf(file = paste0(output_dir, cell_type, "_APOEvsAPOE_AD_upset_", tag, ".pdf"))
    print(upset(
      fromList(sig_genes_ls),
      nsets = length(sig_genes_ls),
      order.by = "freq",
      sets.bar.color = color_map$color,
      text.scale = 1.5
    ))
    dev.off()
  }, error = function(e) {
    message(paste(cell_type, e))
  }
  )
}
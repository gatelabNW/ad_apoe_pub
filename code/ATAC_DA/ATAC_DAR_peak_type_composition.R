# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-11-2022
# Written by: Natalie Piehl
# Summary: Look at peak type composition of DARs
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
})

# Organize inputs
dar_dir <- "/path/to/da/results"
celltype_colors_path <- "/path/to/celltype_color_map.csv"
output_dir <- "/path/to/output/folder"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Define comparison
comparison <- "ADvsHC"

#------------------------------------------------------------------------------
# Make bar plot of composition
#-------------------------------------------------------------------------------

# Get cell types
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2,
                                                function(x) {gsub(" ", "_", x)})
cell_types <- celltype_colors$predicted.celltype.l2

# Create list with sig genes for each cell type
type_ls <- list()
for (cell_type in cell_types) {
  print(cell_type)

  # Load in degs
  tryCatch({
    # Read in peaks
    da_peaks <- read.csv(paste0(dar_dir, comparison, "/", cell_type, "_", comparison, "_dars.csv"))
  
    # Identify sig regions
    da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & abs(da_peaks$avg_log2FC) > lfc.thresh),]
    
    # Calculate proportion of each peak type
    peak_types <- table(da_peaks$peakType)
    
    # Create formatted
    peak_types_formatted <- list(Distal = 0,
                                 Exonic = 0,
                                 Intronic = 0,
                                 Promoter = 0)
    for (type in names(peak_types_formatted)) {
      if (type %in% names(peak_types)) {
        peak_types_formatted[[type]] <- peak_types[[type]]
      }
    }

    # Add proportions to list
    type_ls[[cell_type]] <- as.vector(unlist(peak_types_formatted))
  }, error = function(e) {
    message(e)
  })
}

# Remove any celltypes which don't any dars
type_ls <- type_ls[ lapply(type_ls, sum) > 0 ]

# Convert to dataframe
res <- data.frame(type_ls) %>% t() %>% as.data.frame
colnames(res) <- c("Distal", "Exonic", "Intronic", "Promoter")

# Calculate number of dars
num_dars <- as.vector(rowSums(res))
num_dars_ordered <- num_dars[order(num_dars, decreasing = TRUE)]

# Format
res <- res / rowSums(res)
res$cell_type <- rownames(res)
res$num_dars <- num_dars
res <- res[order(res$num_dars, decreasing = TRUE),]
cell_type_order <- res$cell_type

# Convert to long format
res_long <- res %>%
              pivot_longer(cols = -c(cell_type, num_dars), 
                           names_to = "Peak_Type",
                           values_to = 'Proportion')
res_long$cell_type <- factor(res_long$cell_type, levels = cell_type_order)

# Generate plot
p <- ggplot(res_long, aes(x = cell_type, y = Proportion, fill = Peak_Type)) +
      geom_bar(stat = "identity", color = "black") +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = c("plum1", "violetred1", "skyblue", "palegreen")) +
      theme_Publication_blank() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "right",
            text = element_text(size = 24))
for (i in seq(length(unique(res$cell_type)))) {
  p <- p + annotate("text", x = i, y = 0.95,
                    label = as.character(num_dars_ordered[i]),
                    size = 4.5)
}
p

# Export plot
set_panel_size(p, file = paste0(output_dir, comparison, "_DARs_peak-type_composition.pdf"))
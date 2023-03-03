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
# Written by: Natalie Piehl
# Summary: Make bar plot of peak type composition
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
})

# Organize inputs
celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/broad_celltype_color_map.csv"
da_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/main/out_NP_02-06-2023/"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/composition/out_NP_02-08-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# List comparisons
comparisons <- c("ADvsHC",
                 "ADvsHC_33", "ADvsHC_34", "ADvsHC_44",
                 "34vs33_HC", "34vs33_AD",
                 "44vs33_HC", "44vs33_AD",
                 "44vs34_HC", "44vs34_AD")

#-------------------------------------------------------------------------------
# Run TF enrichment

# Define celltypes
color_map <- read.csv(celltype_colors_path)
cell_types <- sapply(color_map$predicted.celltype.l2,
                     function(x) {gsub(" ", "_", x)}) %>% as.vector

# # Remove platelets cause not enough cells and really skews one of plots
# cell_types <- cell_types[-28]

# Define output list
res_list <- list()

for (comparison in comparisons[c(6,8,10)]) {
  for (cell_type in cell_types) {
    peaks_path <- paste0(da_base_dir, comparison, "/", cell_type, "_", comparison, "_dars.csv")
    if (file.exists(peaks_path)) {
      tryCatch({
        # Read in DARs
        da_peaks <- read.csv(peaks_path)
        
        # Identify up and down dars
        up_dars <- da_peaks[which(da_peaks$BH < padj.thresh & da_peaks$avg_log2FC > lfc.thresh),]
        down_dars <- da_peaks[which(da_peaks$BH < padj.thresh & da_peaks$avg_log2FC < -lfc.thresh),]
        
        # Append to list
        res_list[[cell_type]] <- c(nrow(up_dars), nrow(down_dars))
      }, error = function(e) {
        message(e)
        res_list[[cell_type]] <- c(0, 0)
      })  
    }
  }
  
  # Convert to dataframe
  res <- data.frame(res_list) %>% t() %>% data.frame()
  names(res) <- c("up", "down")
  
  # Create ordering
  res$num_dars <- rowSums(res)
  res <- res[order(res$num_dars, decreasing = TRUE),]
  res <- res[which(res$num_dars != 0),]
  cell_type_order <- rownames(res)
  
  # Calculate proportions
  res_prop <- res[,1:2] / res$num_dars
  
  # Add cell type column
  res_prop$cell_type <- rownames(res_prop)
  
  # Convert to long format
  res_prop_long <- res_prop %>%
    pivot_longer(cols = -c(cell_type), 
                 names_to = "Peak_Type",
                 values_to = 'Proportion')
  res_prop_long$cell_type <- factor(res_prop_long$cell_type, levels = cell_type_order)
  res_prop_long$Peak_Type <- factor(res_prop_long$Peak_Type, levels = c("up", "down"))
  
  # Generate plot
  p <- ggplot(res_prop_long, aes(x = cell_type, y = Proportion, fill = Peak_Type)) +
    geom_bar(stat = "identity", color = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("tomato", "steelblue1")) +
    theme_Publication_blank() +
    # ggtitle(paste0(comparison, " DARs composition")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "right",
          text = element_text(size = 24))
  
  for (i in seq(length(unique(res_prop$cell_type)))) {
    p <- p + annotate("text", x = i, y = 0.95,
                      label = as.character(res$up[i]),
                      size = 4.5)
    p <- p + annotate("text", x = i, y = 0.05,
                      label = as.character(res$down[i]),
                      size = 4.5)
  }
  print(p)
  
  # Export plot
  set_panel_size(p, file = paste0(output_dir, comparison, "_DARs_upvsdown_composition.pdf"))
}

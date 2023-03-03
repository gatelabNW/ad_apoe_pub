# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-19-2023
# Written by: Natalie Piehl
# Summary: Plot TF enrichment in AD vs HC DARs
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
celltype_colors_path <- "/path/to/celltype_color_map.csv"
tf_base_dir <- "/path/to/tf_enrichment/results/"
output_base_dir <- "/path/to/output/folder/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Generate ranked TF enrichment plots
#-------------------------------------------------------------------------------

# Read in colors and celltypes
colors <- read.csv(celltype_colors_path)
cell_types <- gsub(" ", "_", colors$predicted.celltype.l2)

# List comparisons
comparisons <- c("ADvsHC",
                 "ADvsHC_33", "ADvsHC_34", "ADvsHC_44",
                 "34vs33_HC", "34vs33_AD",
                 "44vs33_HC", "44vs33_AD",
                 "44vs34_HC", "44vs34_AD")

for (comparison in comparisons[2:4]) {
  # Define tf dir
  tf_dir <- paste0(tf_base_dir, comparison, "/")
  
  # Generate output dir
  output_dir <- paste0(output_base_dir, comparison, "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (cell_type in cell_types) {
    up_file_name <- paste0(tf_dir, cell_type, "_", comparison, "_enriched_motifs_inUpDARs.csv")
    down_file_name <- paste0(tf_dir, cell_type, "_", comparison, "_enriched_motifs_inDownDARs.csv")
    if (file.exists(up_file_name) & file.exists(down_file_name)) {
      tryCatch({
        # Read in motifs
        up_enriched_motifs <- read.csv(up_file_name)
        down_enriched_motifs <- read.csv(down_file_name)
        
        # Calculate PFC and extract columns of interest
        up_enriched_motifs$LFC <- log2(up_enriched_motifs$fold.enrichment)
        up_enriched_motifs$PFC_up <- abs(up_enriched_motifs$LFC * log10(up_enriched_motifs$p.adjust))
        up_enriched_motifs <- up_enriched_motifs[,c("motif.name", "PFC_up")]
        up_enriched_motifs <- up_enriched_motifs[order(up_enriched_motifs$PFC_up, decreasing = TRUE),]
        up_enriched_motifs$rank_up <- as.numeric(factor(-up_enriched_motifs$PFC_up))
        
        down_enriched_motifs$LFC <- log2(down_enriched_motifs$fold.enrichment)
        down_enriched_motifs$PFC_down <- abs(down_enriched_motifs$LFC * log10(down_enriched_motifs$p.adjust))
        down_enriched_motifs <- down_enriched_motifs[,c("motif.name", "PFC_down")]
        down_enriched_motifs <- down_enriched_motifs[order(down_enriched_motifs$PFC_down, decreasing = TRUE),]
        down_enriched_motifs$rank_down <- as.numeric(factor(-down_enriched_motifs$PFC_down))
        
        # Calculate rank difference
        enriched_motifs <- merge(up_enriched_motifs, down_enriched_motifs)
        enriched_motifs$rank_diff <- enriched_motifs$rank_down - enriched_motifs$rank_up
        write.csv(enriched_motifs, file = paste0(output_dir, cell_type, "_", comparison, "_ADdars_vs_HCdars_enrichment.csv"))

        # Generate plot
        p <- ggplot(enriched_motifs, aes(x = PFC_down, y = PFC_up, fill = rank_diff)) +
          geom_point(size = 2, shape = 21, color = "black") +
          scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
          ggrepel::geom_label_repel(
            data = enriched_motifs, aes(x = PFC_down,
                                        y = PFC_up,
                                        label = motif.name),
            size = 2,
            force = 5,
            color = "black",
            max.overlaps = 15
          ) +
          theme_Publication_blank() +
          ylab("TF enrichment in AD accessible regions (PFC)") +
          xlab("TF enrichment in HC accessible regions (PFC)") +
          ggtitle(paste(cell_type, comparison)) +
          theme(legend.position = "right")

        # Export plot
        set_panel_size(p, file = paste0(output_dir, cell_type, "_", comparison, "_ADdars_vs_HCdars_enrichment.pdf"),
                       width = unit(4, "in"), height = unit(3, "in"))
      }, error = function(e) {
        message(paste(cell_type, e))
      })
    }
  }
}
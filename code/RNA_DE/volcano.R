# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-31-2023
# Written by: Natalie Piehl
# Summary: Make volcano plot with pseudobulk coloring
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("colorspace")
})

# Organize inputs
base_dir <- "/projects/b1169/projects/AD_APOE/results/"
output_base_dir <- "/projects/b1169/projects/AD_APOE/results/de/volcano_with_pseudobulk/out_NP_05-31-2023/"

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

# Define comparison
comparison <- "ADvsHC_33_clonal"
output_dir <- paste0(output_base_dir, comparison, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Define volcano plot function
#-------------------------------------------------------------------------------

# Define volcano plotting function
pseudobulk_volcano_plot <- function(data, pb_degs, file = NULL, title = NULL,
                         padj.lim = NULL, lfc.lim = NULL, lfc.asymmetric = NULL,
                         padj.thresh = 0.01, lfc.thresh = 0.25, x_title = "avg_log2FC",
                         width=unit(4, "inch"), height=unit(4, "inch")) {
  # Create gene name column
  data$gene <- rownames(data)
  
  # Generate PFC scores
  data$PFC <- -log10(data$BH) * abs(data$avg_log2FC)
  PFC <- unique(data$PFC)
  data$PFC[data$PFC == Inf] <- sort(PFC, partial=length(PFC)-1)[length(PFC)-1]
  if (length(PFC) != 1) {
    data$PFC[data$PFC == Inf] <- sort(PFC, partial=length(PFC)-1)[length(PFC)-1]
  }
  
  #Generate log padj column
  data$log_padj <- -log10(data$BH)
  
  # Define limits if not provided
  if (is.null(padj.lim)) {
    log.padj.lim <- unique(data$log_padj)[order(-unique(data$log_padj))][2]
  } else {
    log.padj.lim <- -log10(padj.lim)
  }
  if (is.null(lfc.lim)) {
    lfc.lim <- abs(data[order(-abs(data$avg_log2FC)),"avg_log2FC"][1])
  }
  
  # Generate color column
  data$color <- rep("black", nrow(data))
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC > lfc.thresh), 'color'] <- "lightred"
  data[which(data$BH <= padj.thresh &
               data$avg_log2FC < -lfc.thresh), 'color'] <- "lightblue"
  
  # Add PB coloring
  data[pb_degs, "color"] <- "pb"
  data[which(data$color == "pb" & data$avg_log2FC > lfc.thresh), "color"] <- "red"
  data[which(data$color == "pb" & data$avg_log2FC < -lfc.thresh), "color"] <- "blue"
  
  # Scale down genes outside of bounds
  data[which(data$log_padj > log.padj.lim), 'log_padj'] <- log.padj.lim
  data[which(data$avg_log2FC > lfc.lim), 'avg_log2FC'] <- lfc.lim
  data[which(data$avg_log2FC < -lfc.lim), 'avg_log2FC'] <- -lfc.lim
  
  # Generate asymmetric lfc limits if necessary
  if (is.null(lfc.asymmetric)) {
    lfc_lims <- c(-lfc.lim, lfc.lim)
  } else {
    lfc_lims <- lfc.asymmetric
  }
  
  # Plot data
  p <-
    ggplot(data,
           aes(
             x = avg_log2FC,
             y = log_padj,
             color = color,
             label = gene,
             size = PFC
           )) +
    theme_Publication_blank() +
    geom_hline(
      yintercept = -log10(padj.thresh),
      linetype = 2,
      color = "gray"
    ) +
    geom_vline(xintercept = lfc.thresh,
               linetype = 2,
               color = "gray") +
    geom_vline(
      xintercept = -lfc.thresh,
      linetype = 2,
      color = "gray"
    ) +
    geom_point(aes(size = PFC), alpha = 0.5) +
    scale_color_manual(values = c("red" = "red",
                                  "lightred" = lighten("red", 0.5),
                                  "lightblue" = lighten("blue", 0.5),
                                  "black" = "black",
                                  "blue" = "blue")) +
    geom_text_repel(
      data = data[which(data$color %in% c("red", "blue")), ],
      inherit.aes = T,
      color = 'black',
      size = 4,
      force = 3
    ) +
    theme(legend.position = "none") +
    labs(title = title,
         x = x_title) +
    scale_x_continuous(limits = lfc_lims) +
    scale_y_continuous(limits = c(0, log.padj.lim)) 
  theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12))
  
  # Export plot
  if (is.null(file)) {
    return(p)
  } else {
    set_panel_size(
      p,
      file = file,
      width = width,
      height = height
    )
    return(NULL)
  }
}

#-------------------------------------------------------------------------------
# Generate volcano plot
#-------------------------------------------------------------------------------

# Read in intersections
full_deg_list <- read.csv(paste0(base_dir, "de/MAST_edgeR/out_NP_05-19-2023/full_intersection/", comparison, "_padj0.05_lfc0.125_MAST+edgeR_DEGs.csv"))

# Define dirs
if (!grepl("clonal", comparison) & comparison %!in% c("ADvsHC_44", "44vs34_HC")) {
  de_dir <- paste0(base_dir, "MAST_withEthnicity/out_NP_05-22-2022/", comparison, "/")
} else if (comparison == "ADvsHC_44") {
  de_dir <- paste0(base_dir, "de/diagnosis_44/out_NP_09-28-2022_covarSex/")
} else if (comparison == "44vs34_HC") {
  de_dir <- paste0(base_dir, "de/apoe44vs34_hc/out_NP_11-10-2022_covarSex/")
}

if (grepl("clonal", comparison) & comparison %in% c("ADvsHC_33_clonal", "ADvsHC_34_clonal", "ADvsHC_44_clonal")) {
  de_dir <- paste0("/projects/b1169/projects/AD_APOE/results/tcr/clonality/diagnosisByAPOE_sexEthCovar/out_2023_05_31_AR_NPformatted/", comparison, "/")
} else if (comparison == "44vs34_HC_clonal") {
  de_dir <- paste0("/projects/b1169/projects/AD_APOE/results/tcr/clonality/out_2022_12_14_AR_NPformatted/", comparison, "/")
} else {
  de_dir <- mast_dir <- paste0("/projects/b1169/projects/AD_APOE/results/tcr/clonality/genotype_ethnicity_covar/out_2023_05_30_AR_NPformatted/", comparison, "/")
}

# Define cell types
cell_types <- c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                "CD14_Mono", "CD16_Mono",
                "ASDC", "cDC1", "cDC2", "pDC",
                "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                "NK", "NK_Proliferating", "NK_CD56bright",
                "ILC", "dnT", "gdT")

for (cell_type in cell_types) {
  tryCatch({
    # Read in DEGs
    celltype_files <- list.files(de_dir, pattern = cell_type, full.names = TRUE)
    if (cell_type == "NK") {
      files <- grep(".csv", celltype_files, value = TRUE)
      degs <- read.csv(files[!grepl("CD56|Prolif", files)], row.names = 1)
    } else {
      degs <- read.csv(grep(".csv", celltype_files, value = TRUE), row.names = 1)
    }
    
    # Grab pb degs
    pb_degs <- full_deg_list[,cell_type][!is.na(full_deg_list[,cell_type]) & full_deg_list[,cell_type] != ""]

    # Create volcano plot
    pseudobulk_volcano_plot(degs, pb_degs,
                            title = paste0(comparison, " : ", cell_type),
                            file = paste0(output_dir, cell_type, "_", comparison, "_volcano.pdf"),
                            padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  }, error = function(e) {
    message(e)
  })
}
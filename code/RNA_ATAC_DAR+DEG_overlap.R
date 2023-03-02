# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-10-2022
# Written by: Natalie Piehl
# Summary: Compare DA and DE over AD
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Signac")
  library("Seurat")
})

# Organize inputs
celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/celltype_color_map.csv"
broad_celltype_colors_path <- "/projects/b1169/projects/AD_APOE/data/color/broad_celltype_color_map.csv"
ranges_path <- "/projects/b1169/projects/AD_APOE/data/ranges/full_ranges.rds"
da_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/main/out_NP_02-06-2023/"
output_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_de_comparison/ADvsHC_scatter/out_NP_02-08-2023_maxOver20/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Specify peak type
peak_type <- c("all")

# Specify comparison
comparison <- "ADvsHC_44"

#-------------------------------------------------------------------------------
# Compare DEGs and DARs

# Define de dir
if (comparison == "ADvsHC") {
  de_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis/out_NP_09-28-2022_covarSex+APOE/"
} else if (comparison == "ADvsHC_33") {
  de_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis_33/out_NP_09-28-2022_covarSex/"
} else if (comparison == "ADvsHC_34") {
  de_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis_34/out_NP_09-28-2022_covarSex/"
} else if (comparison == "ADvsHC_44") {
  de_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis_44/out_NP_09-28-2022_covarSex/"
}

# Define specific output dir
output_dir <- paste0(output_base_dir, comparison, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read in ranges
ranges <- readRDS(ranges_path)

# Get fine cell types
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors$predicted.celltype.l2 <- sapply(celltype_colors$predicted.celltype.l2,
                                                function(x) {gsub(" ", "_", x)})
cell_types <- celltype_colors$predicted.celltype.l2

# Get broad cell types
broad_celltype_colors <- read.csv(broad_celltype_colors_path)
broad_celltype_colors$predicted.celltype.l2 <- sapply(broad_celltype_colors$predicted.celltype.l2,
                                                function(x) {gsub(" ", "_", x)})
broad_cell_types <- broad_celltype_colors$predicted.celltype.l2

# Make broad cell type map
broad_map <- data.frame(fine = c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                                "CD14_Mono", "CD16_Mono",
                                "ASDC", "cDC1", "cDC2", "pDC",
                                "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                                "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                                "NK", "NK_Proliferating", "NK_CD56bright",
                                "ILC", "dnT", "gdT",
                                "Platelet", "Eryth", "HSPC", "Doublet"),
                         broad = c(rep("B_Cells", 4),
                                rep("Monocytes", 2),
                                rep("Dendritic_Cells", 4),
                                rep("CD4+_T_Cells", 6),
                                rep("CD8+_T_Cells", 5),
                                rep("NK_Cells", 3),
                                rep("Other_T_Cells", 3),
                                rep("Other", 4)))

for (cell_type in cell_types) {
  tryCatch({
    # Load DEGs
    degs <- read.csv(paste0(de_dir, cell_type, "_AD_vs_HC_degs.csv"))
    degs <- degs[,c(1,3,7)]
    names(degs) <- c("nearestGene", "RNA_LFC", "RNA_BH")
    
    # Load DARs
    broad_cell_type <- broad_map[which(broad_map$fine == cell_type), "broad"]
    dars <- read.csv(paste0(da_dir, comparison, "/", broad_cell_type, "_", comparison, "_dars.csv"))
    if ("all" %!in% peak_type) {
      dars <- dars[which(dars$peakType %in% peak_type),]
    }
    dars <- dars[,c(2,4,8,9,10,11)]
    names(dars)[2:3] <- c("ATAC_LFC", "ATAC_BH")
    
    # Merge together
    res <- merge(degs, dars)
    
    # Add color column
    res$color <- rep("notSig", nrow(res))
    res[which((abs(res$RNA_LFC) > lfc.thresh & res$RNA_BH < padj.thresh) | (abs(res$ATAC_LFC) > lfc.thresh & res$ATAC_BH < padj.thresh)), "color"] <- "sigInOnlyOne"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "UpATAC+UpRNA"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "DownATAC+DownRNA"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "DownATAC+UpRNA"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "UpATAC+DownRNA"
    
    # Calculate limits for square visual
    lim <- max(max(abs(res$RNA_LFC)), max(abs(res$ATAC_LFC)))
    
    # Plot data
    p <- ggplot(res,
             aes(
               x = ATAC_LFC,
               y = RNA_LFC,
               fill = color,
               label = gene_region_index,
               size = abs(ATAC_LFC*RNA_LFC)
             )) +
      theme_Publication_blank() +
      theme_minimal() +
      geom_point(alpha = 0.5, shape = 21, stroke = 0.2) +
      geom_vline(aes(xintercept = 0)) +
      geom_hline(aes(yintercept = 0)) +
      geom_vline(xintercept = lfc.thresh,
                 linetype = 2,
                 color = "gray") +
      geom_vline(xintercept = -lfc.thresh,
                 linetype = 2,
                 color = "gray") +
      geom_hline(yintercept = lfc.thresh,
                 linetype = 2,
                 color = "gray") +
      geom_hline(yintercept = -lfc.thresh,
                 linetype = 2,
                 color = "gray") +
      scale_fill_manual(values = c("notSig" = "gray90",
                                   "sigInOnlyOne" = "gray30",
                                    "UpATAC+UpRNA" = "red",
                                    "DownATAC+DownRNA" = "blue",
                                    "UpATAC+DownRNA"= "green",
                                    "DownATAC+UpRNA"= "green")) +
      geom_label_repel(
        data = res[which(res$color %in% c("UpATAC+UpRNA", "DownATAC+DownRNA",
                                          "UpATAC+DownRNA", "DownATAC+UpRNA")), ],
        inherit.aes = T,
        fill = 'white',
        color = 'black',
        max.overlaps = 20,
        size = 2,
        force = 5
      ) +
      labs(title = paste0(cell_type, " DE vs ", broad_cell_type, " DA ", comparison)) +
      scale_size_continuous(range = c(0.001, 3)) +
      scale_x_continuous(limits = c(-lim, lim)) +
      scale_y_continuous(limits = c(-lim, lim)) +
      theme(axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(colour = NA, fill = "transparent"))
    print(p)
    
    # Export plot and results
    set_panel_size(p, file = paste0(output_dir, cell_type, "-DEvs_", broad_cell_type, "-DA_LFCscatter.pdf"),
                   width = unit(3, "in"), height = unit(3, "in"))
    write.csv(res, paste0(output_dir, cell_type, "-DEvs_", broad_cell_type, "-DA.csv"))
    print(paste(cell_type, "completed!"))
  }, error = function(e) {
    message(e)
  })
}
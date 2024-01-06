# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-03-2022
# Written by: Natalie Piehl
# Summary: Make heatmap of celltype markers
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
})

# Initialize input parameters
celltype_colors_path <- "/path/to/celltype_color_map.csv"
seurat_object <- "/path/to/object"
output_dir <- "/path/to/directory/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Data formatting

# Load Seurat object
load(seurat_object)

# Get colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[seq(24,nrow(celltype_colors)+2),] <- celltype_colors[seq(22,nrow(celltype_colors)),]
celltype_colors[22:23,] <- celltype_colors[28:29,]
celltype_colors <- celltype_colors[-c(28:29),]

# Order colors in seurat object
celltype_colors <- celltype_colors[which(celltype_colors$predicted.celltype.l2 %in% unique(s[["predicted.celltype.l2"]])[,1]),]
s@meta.data$predicted.celltype.l2 <- factor(s@meta.data$predicted.celltype.l2,
                                            levels = celltype_colors$predicted.celltype.l2)

# Cell type markers
markers<-c("CD3D","CD4","CD8A","CD8B","NCAM1","KLRB1","NKG7","CD19",
           "MS4A1","CD38","CD14","CD27","CD68","FCGR3A","CD1C","LILRA4",
           "PPBP","CD34","JCHAIN","HBB",
           # "CD44","SELL",
           "MKI67","FOXP3","IL7R","CCL4L2",
           "S100A8", "S100A9", "GNLY", "GZMB", "HLA-DRA", "TCF7","LEF1","CCR7")

# Sample maximum cells from each celltype
barcode_celltype <- s[["predicted.celltype.l2"]]
barcode_celltype$barcode <- rownames(barcode_celltype)
barcode_sampled <- barcode_celltype %>%
  group_by(predicted.celltype.l2) %>%
  slice_sample(n = 500, replace = TRUE)
s_sub <- subset(s, cells = unique(barcode_sampled$barcode))

# Visualize heatmap of cell type markers
s <- SetIdent(s_sub, value = "predicted.celltype.l2")
myplot <- DoHeatmap(s,
                    group.colors = celltype_colors$new_color,
                    features = markers,
                    slot = "data",
                    assay = "RNA"
) +
  scale_fill_gradientn(colors = c("white", "darkorchid2", "navy"))

set_panel_size(myplot,
               file = paste0(output_dir, "rna_celltype_markers_heatmap.pdf"),
               width = unit(12, "in"), height = unit(6, "in")
)

# Visualize with dotplot
s <- SetIdent(s, value = "predicted.celltype.l2")
myplot <- DotPlot(s,
                  cols = c("white", "darkorchid2", "navy"),
                  features = markers
) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_discrete(limits=rev)

set_panel_size(myplot,
               file = paste0(output_dir, "rna_celltype_markers_dotplot.pdf"),
               width = unit(8, "in"), height = unit(8, "in")
)
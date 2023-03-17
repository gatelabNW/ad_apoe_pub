# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-30-2022
# Written by: Natalie Piehl
# Summary: Identify cell types using CITEseq reference based mapping
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("SeuratDisk")
})

# Initialize input parameters
seurat_object <- "/path/to/integrated/seurat/object"
ref_atlast_path <- "/path/to/pbmc_CITEseq_ref"
output_dir <- "/path/to/output/folder"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Run supervised clustering

# Load in seurat object
load(seurat_object)
s

# Load in reference
reference <- LoadH5Seurat(ref_atlast_path)

# Set assay to SCT
DefaultAssay(object = s) <- "SCT"

# Find transfer anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = s,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

# Map to reference atlas
s <- MapQuery(
  anchorset = anchors,
  query = s,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)
s

# Save updated object
save(s, file = paste0(output_dir, "s_sup_clustering"))

# Order cell type colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors <- celltype_colors[which(celltype_colors$predicted.celltype.l2 %in% unique(s[["predicted.celltype.l2"]])[,1]),]
s$predicted.celltype.l2 <- factor(s$predicted.celltype.l2,
                                  levels = celltype_colors$predicted.celltype.l2)

# Generate original UMAP plots
p1 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
p2 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE,
             cols = celltype_colors$color)
set_panel_size(p1, file = paste0(output_dir, "og_umap_L1.pdf"), width = unit(5, "in"), height = unit(4, "in"))
set_panel_size(p2, file = paste0(output_dir, "og_umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))

# Generate ref UMAP plot
p2 = DimPlot(s, reduction = "ref.umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE,
             cols = celltype_colors$color)
set_panel_size(p2, file = paste0(output_dir, "umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))

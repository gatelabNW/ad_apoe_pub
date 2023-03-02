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
# Summary: Run supervised clustering
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
seurat_object <- "/projects/b1169/projects/AD_APOE/results/seurat/clustering_no1086/out_NP_08-24-2022/s_integrated"
# seurat_object <- "/projects/b1169/nat/CSF_immunity_age/results/tcr_cleanup/main/out/s_tcrclean"
ref_atlast_path <- "/projects/b1169/nat/CSF_immunity_age/data/celltype_atlas/pbmc_multimodal.h5seurat"
output_dir <- "/projects/b1169/projects/AD_APOE/results/seurat/supervised_clustering/out_NP_08-30-2022/"
# test_seurat_object <- paste0(output_dir, "s_test.rds")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#------------------------------------------------------------------------------
# Run supervised clustering

# Load in seurat object
load(seurat_object)
s

# Idents(s) <- "orig.ident"
# samples_to_use <- unique(s[["orig.ident"]])[1:10,1]
# s <- subset(s, orig.ident %in% samples_to_use)
# s <- subset(s, downsample = 500)

# Load in reference
reference <- LoadH5Seurat(ref_atlast_path)

# Run standard normalization
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

# Visualize results
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"
celltype_colors <- celltype_colors[which(celltype_colors$predicted.celltype.l2 %in% unique(s[["predicted.celltype.l2"]])[,1]),]
s$predicted.celltype.l2 <- factor(s$predicted.celltype.l2,
                                  levels = celltype_colors$predicted.celltype.l2)
p1 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
p2 = DimPlot(s, reduction = "umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE,
             cols = celltype_colors$color)
set_panel_size(p1, file = paste0(output_dir, "og_umap_L1.pdf"), width = unit(5, "in"), height = unit(4, "in"))
set_panel_size(p2, file = paste0(output_dir, "og_umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))



p2 = DimPlot(s, reduction = "ref.umap", group.by = "predicted.celltype.l2",
             label = TRUE, label.size = 3 ,repel = TRUE,
             cols = celltype_colors$color)
set_panel_size(p2, file = paste0(output_dir, "umap_L2.pdf"), width = unit(5, "in"), height = unit(4, "in"))
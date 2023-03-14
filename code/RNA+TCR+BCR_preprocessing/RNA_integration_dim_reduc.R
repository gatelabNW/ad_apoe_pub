# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-18-2022
# Written by: Natalie Piehl
# Summary: Integrate and cluster RNA data
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
seurat_object <- "/path/to/RNA/seurat/object"
output_dir <- "/path/to/output/folder/"
pca_dir <- paste0(output_dir, "pca/")
umap_dir <- paste0(output_dir, "umap/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(umap_dir, showWarnings = FALSE, recursive = TRUE)

# Set parameters
pc_dim <- 15
res_thresh <- 0.3
num_features <- 1000
n_neighbors <- 20
min_dist <- 0.4
n_epochs <- 250

#------------------------------------------------------------------------------
# Split, transform, then reintegrate data sets (SCT normalization)

# Load Seurat object
load(seurat_object)

# Remove 1086 and 659
s
s <- subset(s, orig.ident %!in% c("G1086_y4", "G659_y2"))
s

# Split by ID
s_split <- SplitObject(s, split.by = "orig.ident")

# Normalize each sample separately
for (i in names(s_split)) {
  s_split[[i]] <- SCTransform(s_split[[i]],
                              variable.features.n = num_features,
                              ncells = 10000,
                              vars.to.regress = c(
                                "nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt"
                              )
  )
  DefaultAssay(s_split[[i]]) <- "SCT"
}

# Save SCTransformed objects
save(s_split, file = paste0(output_dir, "s_split_sctransformed"))
# load(paste0(output_dir, "s_split_sctransformed"))

# Prep for integration
integration_features <- SelectIntegrationFeatures(object.list = s_split,
                                                  assay = rep("SCT", length(s_split)))
s_split <- PrepSCTIntegration(object.list = s_split,
                              anchor.features = integration_features)

# Select references
reference_dataset <- which(names(s_split) %in% c("G1162", "G1055_y2"))

# Find anchors
anchors <- FindIntegrationAnchors(object.list = s_split,
                                  normalization.method = "SCT",
                                  anchor.features = integration_features,
                                  reference = reference_dataset
                                  )

# Integrate data
s <- IntegrateData(anchorset = anchors,
                   normalization.method = "SCT")
s

# Checkpoint save/load
save(s, file = paste0(output_dir, "s_integrated"))
# load(paste0(output_dir, "s_integrated"))

#------------------------------------------------------------------------------
# Generate PCA

# Set assay to integrated
DefaultAssay(object = s) <- "integrated"

# Generate PCA
s <- RunPCA(object = s)

# Initialize identities to look at on PCA
idents <- c("orig.ident", "Age", "Sex", "Sort_day", "Diagnosis")

# Plot PCA for each identity
for (ident in idents) {
  s <- SetIdent(s, value = ident)
  myplot <- DimPlot(s, reduction = "pca", pt.size = 2, shuffle = TRUE)
  set_panel_size(myplot,
                 file = paste0(pca_dir, "PCA_", ident, ".pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
}

# Visualize PC Dimensionality of Dataset
myplot <- ElbowPlot(s)
set_panel_size(myplot,
               file = paste0(output_dir, "Elbowplot.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

#------------------------------------------------------------------------------
# UMAP visualization

# Cluster Cells on SCT integrated assay
s <- FindNeighbors(s, dims = 1:pc_dim)
s <- FindClusters(s, resolution = res_thresh)

# Run UMAP on SCT integrated assay
s <- RunUMAP(s,
             dims = 1:pc_dim,
             n.neighbors = n_neighbors,
             min.dist = min_dist,
             n.epochs = n_epochs)

# List identities
idents <- c(idents, "seurat_clusters")

# Visualize each identity on UMAP
for (ident in idents) {
  s <- SetIdent(s, value = ident)
  myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1, shuffle = TRUE)
  set_panel_size(myplot,
                 file = paste0(umap_dir, "UMAP_", ident, ".pdf"),
                 width = unit(5, "in"), height = unit(4, "in")
  )
}

# View Diagnosis split
s <- SetIdent(s, value = "Diagnosis")
myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1,
                  shuffle = TRUE, split.by = "Diagnosis")
set_panel_size(myplot,
               file = paste0(umap__dir, "UMAP_Diagnosis_split.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

# View Sex split
s <- SetIdent(s, value = "Sex")
myplot <- DimPlot(s, reduction = "umap", label = TRUE, pt.size = 1,
                  shuffle = TRUE, split.by = "Sex")
set_panel_size(myplot,
               file = paste0(umap__dir, "UMAP_Sex_split.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)
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
# Written by: Abhi Ramakrishnan
# Summary: Run ArchR (dev branch) pipeline on scATAC data
# Create UMAP and annotate cell types for scATAC data
#
#-------------------------------------------------------------------------------
# Initialization 
# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ArchR")
  library("doMC")
  library("qpdf")
  library("harmony")
  library("rlang")
  library("devtools")
  library("DESeq2")
  library("presto")
  library("HDF5Array")
  library("hexbin")
  library("rhandsontable")
  library("ComplexHeatmap")
})

# Set random seed
set.seed(123)

# Set number of threads
addArchRThreads(threads = 26)

# Define project directory
proj_dir <- "path/to/project/directory/"

# Generate directory for manual outputs outside of project folder
output_dir <- "path/to/output/directory/"

ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)

# Set working dir to project dir
setwd(proj_dir)

# Load project
proj <- loadArchRProject(proj_dir)

# Save project in new directory before next step
# DO NOT load in this copy
proj2_dir <- "path/to/preprocessed/project/"
proj2 <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj2_dir,
                         load = FALSE)
rm(proj2)
#-------------------------------------------------------------------------------
# 1) Dimension reduction
#-------------------------------------------------------------------------------

# Compute LSI (reduced the resolution and increased by 2 iterations compared to default)
proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI",
                        iterations = 4,
                        # force = TRUE, # comment this out and change name if you want to make a new version
                        clusterParams = list(
                          resolution = c(0.2, 0.3, 0.5),
                          sampleCells = 10000,
                          n.start= 10)
)

# Identify clusters (Louvain)
proj <- addClusters(input = proj,
                    name = "Clusters25",
                    reducedDims = "IterativeLSI",
                    maxClusters = 25)

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
#-------------------------------------------------------------------------------
# 2) Generate UMAPs
# Final umap used for further analysis was "UMAP_n20_d0.1"
#-----------------------------------------------------------------------------
# 1.) UMAP with n_neighbors = 20 and Min_dist = 0.1
proj <- addUMAP(ArchRProj = proj,
                reducedDims = "IterativeLSI",
                name = "UMAP_n20_d0.1",
                nNeighbors = 20,
                minDist = 0.1, )
# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
# Generate plots
p <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Clusters25",
                   embedding = "UMAP_n20_d0.1")
q <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Sample",
                   embedding = "UMAP_n20_d0.1")


# Export plots
plotPDF(p,
        name = "umap_n20_d0.1_clusters25.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
plotPDF(q,
        name = "umap_n20_d0.1_sample.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)


# 2.) UMAP with n_neighbors = 20 and Min_dist = 0.2
proj <- addUMAP(ArchRProj = proj,
                reducedDims = "IterativeLSI",
                name = "UMAP_n20_d0.2",
                nNeighbors = 20,
                minDist = 0.2, )


# Generate plots
p <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Clusters25",
                   embedding = "UMAP_n20_d0.2")
q <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Sample",
                   embedding = "UMAP_n20_d0.2")
r <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "predictedGroup",
                   embedding = "UMAP_n20_d0.2")

# Export plots
plotPDF(p,
        name = "umap_n20_d0.2_clusters25.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5
        )
plotPDF(q,
        name = "umap_n20_d0.2_sample.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
plotPDF(r,
        name = "umap_n20_d0.2_predictedGroup.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)

# 3.) UMAP with n_neighbors = 30 and Min_dist = 0.3
proj <- addUMAP(ArchRProj = proj,
                reducedDims = "IterativeLSI",
                name = "UMAP_n30_d0.3",
                nNeighbors = 30,
                minDist = 0.3, )


# Generate plots
p <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Clusters25",
                   embedding = "UMAP_n30_d0.3")
q <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Sample",
                   embedding = "UMAP_n30_d0.3")
r <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "predictedGroup",
                   embedding = "UMAP_n30_d0.3")

# Export plots
plotPDF(p,
        name = "umap_n30_d0.3_clusters25.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5
)
plotPDF(q,
        name = "umap_n30_d0.3_sample.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
plotPDF(r,
        name = "umap_n30_d0.3_predictedGroup.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
# 4.) UMAP with n_neighbors = 30 and Min_dist = 0.2
proj <- addUMAP(ArchRProj = proj,
                reducedDims = "IterativeLSI",
                name = "UMAP_n30_d0.2",
                nNeighbors = 30,
                minDist = 0.2, )


# Generate plots
p <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Clusters25",
                   embedding = "UMAP_n30_d0.2")
q <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Sample",
                   embedding = "UMAP_n30_d0.2")
r <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "predictedGroup",
                   embedding = "UMAP_n30_d0.2")

# Export plots
plotPDF(p,
        name = "umap_n30_d0.2_clusters25.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5
)
plotPDF(q,
        name = "umap_n30_d0.2_sample.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
plotPDF(r,
        name = "umap_n30_d0.2_predictedGroup.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5)
# Add Harmony to reduce batch effect
proj<- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)

# Generate umaps grouped by sample metadata
v <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Age",
                   embedding = "UMAP_n20_d0.1",
                   continuousSet = "coolwarm")
w <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Diagnosis",
                   embedding = "UMAP_n20_d0.1")
x <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "APOE_Genotype",
                   embedding = "UMAP_n20_d0.1")
y <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "Sex",
                   embedding = "UMAP_n20_d0.1")
z <- plotEmbedding(ArchRProj = proj,
                   colorBy = "cellColData",
                   name = "RNA",
                   embedding = "UMAP_n20_d0.1")

# Manual UMAP cell type annotation using gene score matrix
#-------------------------------------------------------------------------------

markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters25",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveHDF5SummarizedExperiment(markersGS, dir=output_dir, prefix="markersGS", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)

# Plot heatmaps of PBMC cell types marker genes

markerGenes <-c("CD3D","CD4","CD8A","CD8B","NCAM1","KLRB1","NKG7","CD19",
           "MS4A1","CD38","CD14","CD27","CD68","FCGR3A","CD1C","LILRA4",
           "PPBP","CD34","JCHAIN","HBB",
           # "CD44","SELL",
           "MKI67","FOXP3","IL7R","CCL4L2",
           "S100A8", "S100A9", "GNLY", "GZMB", "HLA-DRA", "TCF7","LEF1","CCR7")

gene_list <- getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
)
# Find additional markers in gene score matrix to add to markerGenes
filteredMarkerGenes <- c()
for (gene in markerGenes) {
  if (gene %in% gene_list) {
    filteredMarkerGenes <- c(filteredMarkerGenes, gene)
  }
}
print(filteredMarkerGenes)

# Plot gene score heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.001 & Log2FC >= 0.5",
  labelMarkers = filteredMarkerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap",
        width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

# Use MAGIC to impute gene scores (overcome sparsity of ATAC data)
proj <- addImputeWeights(proj,
                         reducedDims = "Harmony")

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
v <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP_n20_d0.1",
  imputeWeights = getImputeWeights(proj)
)
plotPDF(plotList = v,
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)

#Rearrange for grid plotting
v2 <- lapply(v, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
grid_v2 <- do.call(cowplot::plot_grid, c(list(ncol = 3),v2))
plotPDF(grid_v2,
        name = "Grid-UMAP-Marker-Genes-W-Imputation.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 8, height = 8)
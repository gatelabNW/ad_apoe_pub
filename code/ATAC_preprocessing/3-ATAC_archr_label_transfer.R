# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-22-2022
# Written by: Abhi Ramakrishnan
# Summary: Run ArchR (dev branch) pipeline on scATAC data
# Transfer labels to scATAC-seq data from matched scRNA-seq data
#
#-------------------------------------------------------------------------------
# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ArchR")
  library("doMC")
  library("qpdf")
  library("harmony")
  library("devtools")
  library("DESeq2")
  library("presto")
  library("HDF5Array")
  library("hexbin")
  library("rhandsontable")
  library("Seurat")
  library("ggrastr")
})

# Set random seed
set.seed(123)

# Set number of threads
addArchRThreads(threads = 1)

# Source helper functions
source("path/to/helper_functions.R")

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


# Save project in new directory
# DO NOT load in this copy
proj3_dir <- "path/to/umapped/project/"
proj3 <- saveArchRProject(ArchRProj = proj,
                          outputDirectory = proj3_dir,
                          load = FALSE)
rm(proj3)

#-------------------------------------------------------------------------------
# 1) Unconstrained integration using matched scRNA-seq data
# Perform this step alone and assess results before defining parameters for
# constrained integration
#-------------------------------------------------------------------------------
# Load in scRNA dataset
scRNA <- "path/to/seurat/object/"
load(scRNA)
DefaultAssay(s) <- "RNA"
# Unconstrained integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Un",
  reducedDims = "IterativeLSI",
  seRNA = s,
  addToArrow = FALSE,
  groupRNA = "predicted.celltype.l2",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
)

# Find out which scRNA cell type is most abundant in each scATAC cluster

cM <- as.matrix(confusionMatrix(proj$Clusters25, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
matched_clusters <- cbind(preClust, rownames(cM)) #Assignments
print(matched_clusters)
write.csv(matched_clusters, file= paste0(output_dir, "unconstrained_integration_25_clusters.csv"))

# Visualize scATAC umap using unconstrained integration labels
pal <- paletteDiscrete(values = s$predicted.celltype.l2)
p1 <- plotEmbedding(
  proj,
  colorBy = "cellColData",
  name = "predictedGroup_Un",
  pal = pal,
  embedding = "UMAP_n20_d0.1"
)
plotPDF(p1,
        name = "umap_n20_d0.1_predictedGroup_Un.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5
)

# Save results of unconstrained integration as a separate object
tmp <- proj@cellColData$predictedGroup_Un
save(tmp, file = paste0(output_dir, "annotated_celltypes_Un"))
# Make a table of which cells are in each ATAC cluster and predicted group
clusters25 <- proj$Clusters25
predictedGroup_Un_x_Clusters25 <- table(tmp, clusters25)
write.csv(predictedGroup_Un_x_Clusters25, file = paste0(output_dir, "predictedGroup_Un_x_Clusters25.csv"))

# ---------------------------------------------------------------------
# 2) Constrained Integration
# ---------------------------------------------------------------------
# Save unconstrained integration project in new directory
# DO NOT load in this copy
proj4_dir <- "path/to/unconstrained_integration/project/"
proj4 <- saveArchRProject(ArchRProj = proj,
                          outputDirectory = proj4_dir,
                          load = FALSE)
rm(proj4)
warnings()

# The following constraints use our knowledge of known cell type markers and gene score matrix values for ATAC data
# Make 4 groups of ATAC clusters which are likely to correspond to certain cell types:
# Note: These lists do not have any overlap in ATAC clusters to avoid assigning an ATAC cell to more than one RNA cell

# Load in scRNA dataset
scRNA <- "path/to/seurat/object/"
load(scRNA)
DefaultAssay(s) <- "RNA"
# Make a list of likely T cell and NK cell clusters in scATAC data
clustTNK <- c("C11",  "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19",  "C20", "C21")

# Clusters with B and T cells
clustBT <- c("C7", "C9", "C10")

# Make a list of likely monocyte, B, DC, and other cell clusters in scATAC data
clustBMono <- c("C1", "C2", "C3", "C4", "C5","C6","C8")

# Clusters with a mix of cells
clustAll <- c("C22", "C23", "C24", "C25")

# Similarly, make matching lists of all cell groups in RNA data
cTNK <- paste0(c("CD4 CTL",
                 "CD4 Naive",
                 "CD4 Proliferating",
                 "CD4 TCM",
                 "CD4 TEM",
                 "Treg",
                 "CD8 Naive",
                 "CD8 Proliferating",
                 "CD8 TCM",
                 "CD8 TEM",
                 "NK",
                 "NK Proliferating",
                 "NK_CD56bright",
                 "MAIT",
                 "dnT",
                 "ILC",
                 "gdT"), collapse="|")
cBT <- paste0(c("B intermediate",
                "B memory",
                "B naive",
                "Plasmablast",
                "Eryth",
                "CD4 CTL",
                "CD4 Naive",
                "CD4 Proliferating",
                "CD4 TCM",
                "CD4 TEM",
                "Treg",
                "CD8 Naive",
                "CD8 Proliferating",
                "CD8 TCM",
                "CD8 TEM",
                "NK",
                "NK Proliferating",
                "NK_CD56bright",
                "MAIT",
                "dnT",
                "ILC",
                "gdT"), collapse="|")

cBMono <- paste0(c("ASDC",
                   "cDC2",
                   "pDC",
                   "CD14 Mono",
                   "CD16 Mono",
                   "HSPC",
                   "Platelet",
                   "B intermediate",
                   "B memory",
                   "B naive",
                   "Plasmablast"), collapse="|")

cAll <- paste0(c("CD4 CTL",
                 "CD4 Naive",
                 "CD4 Proliferating",
                 "CD4 TCM",
                 "CD4 TEM",
                 "Treg",
                 "CD8 Naive",
                 "CD8 Proliferating",
                 "CD8 TCM",
                 "CD8 TEM",
                 "NK",
                 "NK Proliferating",
                 "NK_CD56bright",
                 "MAIT",
                 "dnT",
                 "ILC",
                 "gdT",
                 "B intermediate",
                 "B memory",
                 "B naive",
                 "Plasmablast",
                 "Eryth",
                 "ASDC",
                 "cDC2",
                 "pDC",
                 "CD14 Mono",
                 "CD16 Mono",
                 "HSPC",
                 "Platelet"
                 ), collapse="|")

rnaTNK <- colnames(s)[grep(cTNK, s$predicted.celltype.l2)]
head(rnaTNK)
rnaBT <- colnames(s)[grep(cBT, s$predicted.celltype.l2)]
head(rnaBT)
rnaBMono <- colnames(s)[grep(cBMono, s$predicted.celltype.l2)]
head(rnaBMono)
rnaAll <- colnames(s)[grep(cAll, s$predicted.celltype.l2)]
head(rnaAll)

# Make a nested list that identifies cells in each group across both platforms
groupList <- SimpleList(
  TNK = SimpleList(
    ATAC = proj$cellNames[proj$Clusters25 %in% clustTNK],
    RNA = rnaTNK
  ),
  BT = SimpleList(
    ATAC = proj$cellNames[proj$Clusters25 %in% clustBT],
    RNA = rnaBT
  ),
  BMono = SimpleList(
    ATAC = proj$cellNames[proj$Clusters25 %in% clustBMono],
    RNA = rnaBMono
  ),
  All = SimpleList(
    ATAC = proj$cellNames[proj$Clusters25 %in% clustAll],
    RNA = rnaAll
  )
)

# FINALIZED Constrained Integration (4 groups used)
# Pass this group list to groupList parameter to constrain the integration
# Save to arrow files this time
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = s,
  addToArrow = TRUE,
  groupList = groupList,
  groupRNA = "predicted.celltype.l2",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE # added this since values did not save to arrow files
)
# Find out which scRNA cell type is most abundant in each scATAC cluster
proj$predictedGroup
cM2 <- as.matrix(confusionMatrix(proj$Clusters25, proj$predictedGroup))
preClust2 <- colnames(cM2)[apply(cM2, 1 , which.max)]
matched_clusters2 <- cbind(preClust2, rownames(cM2)) #Assignments
print(matched_clusters2)
write.csv(matched_clusters2, file="constrained_integration_final.csv")
# Save results of unconstrained integration as a separate object
tmp <- proj@cellColData$predictedGroup
save(tmp, file = paste0(output_dir, "annotated_celltypes_final"))
# Make a table of which cells are in each ATAC cluster and predicted group
clusters25 <- proj$Clusters25
predictedGroup_x_Clusters25 <- table(tmp, clusters25)
write.csv(predictedGroup_x_Clusters25, file = paste0(output_dir, "predictedGroup_x_Clusters25.csv"))
# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
# See if gene integration matrix is added
getAvailableMatrices(proj)
table(proj$predictedGroup)
print("Gene Integration finished!")

# ---------------------------------------------------------------------
# Visualize scATAC umap using constrained integration labels
# ---------------------------------------------------------------------
pal <- paletteDiscrete(values = s$predicted.celltype.l2)
# Visualize umaps using predictedGroup labels vs ClustersRNA labels

p1 <- plotEmbedding(
  proj,
  colorBy = "cellColData",
  name = "predictedGroup",
  pal = pal,
  embedding = "UMAP_n20_d0.1"
)
plotPDF(p1,
        name = "umap_n20_d0.1_predictedGroup.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 5, height = 5
)
#-------------------------------------------------------------------------------
# Plot Seurat style umap of ATAC data with new RNA cell type labels (modified Natalie's conversion-archr_to_seurat script)
#-------------------------------------------------------------------------------
# Extract UMAP coords
umap_coords <- proj@embeddings$UMAP_n20_d0.1$df
labels <- proj$predictedGroup
# Merge coords and metadata
data <- cbind(umap_coords, labels)
# Load in color map
colormap_path <- "/path/to/celltype/colors.csv"
color_map <- read.csv(colormap_path)
color_map <- color_map[which(color_map$predicted.celltype.l2 %in% data$labels),]
data$labels <- factor(data$labels, levels = color_map$predicted.celltype.l2)
# Generate plot
p <- ggplot(data, aes(x = `IterativeLSI#UMAP_Dimension_1`,
                      y = `IterativeLSI#UMAP_Dimension_2`,
                      color = labels)) +
  scale_color_manual(values = color_map$new_color) +
  theme_Publication_blank() +
  geom_point(size = 0.1) +
  guides(color = guide_legend(override.aes = list(size=5)))
p

# Export plot
set_panel_size(rasterise(p, dpi = 600), file = paste0(output_dir, "atac_UMAP_celltype.pdf"),
               width = unit(6, "in"), height = unit(6, "in"))

# Try it for umap_n30_d0.2
# Extract UMAP coords
umap_coords <- proj@embeddings$UMAP_n30_d0.2$df
labels <- proj$predictedGroup
# Merge coords and metadata
data <- cbind(umap_coords, labels)
# Load in color map
colormap_path <- "/path/to/celltype/colors.csv"
color_map <- read.csv(colormap_path)
color_map <- color_map[which(color_map$predicted.celltype.l2 %in% data$labels),]
data$labels <- factor(data$labels, levels = color_map$predicted.celltype.l2)
# Generate plot
p <- ggplot(data, aes(x = `IterativeLSI#UMAP_Dimension_1`,
                      y = `IterativeLSI#UMAP_Dimension_2`,
                      color = labels)) +
  scale_color_manual(values = color_map$new_color) +
  theme_Publication_blank() +
  geom_point(size = 0.1) +
  guides(color = guide_legend(override.aes = list(size=5)))
p

# Export plot
set_panel_size(rasterise(p, dpi = 600), file = paste0(output_dir, "atac_UMAP_celltype_n30_d0.2.pdf"),
               width = unit(6, "in"), height = unit(6, "in"))

#-------------------------------------------------------------------------------
# Manual UMAP cell type annotation using gene integration matrix
#-------------------------------------------------------------------------------
# Select only cell type groups with enough cells (remove ASDC and HSPC with 3 cells each)
predicted_groups <- data.frame(table(proj$predictedGroup))
selected_groups <- levels(predicted_groups$Var1)
selected_groups <- selected_groups[-c(1,19)]
markersGI <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
  groupBy = "predictedGroup",
  useGroups = selected_groups,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveHDF5SummarizedExperiment(markersGI, dir=output_dir, prefix="markersGI", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)

# Plot heatmaps of cluster marker genes
markerGenes  <- c(
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A", "CD4", "IL7R", #TCells
  "NKG7", #NK Cells
  "CD68" # Dendritic cells
  )
dc_genes <- c(
  "CCR7", "CD45RA", "CD209", "CLEC4C", "NRP1","B220", # pDCs
  "CD163", "CLEC10A", "NOTCH2", "ITGAM", "SIRPA", "CX3CR1", "CD1C", "CD2",  "FCER1A", "CD11C" # cDC2s
)
gene_list <- getFeatures(
  ArchRProj = proj,
  useMatrix = "GeneIntegrationMatrix",
)
# Find additional markers in gene score matrix to add to markerGenes
for (gene in dc_genes) {
  if (gene %in% gene_list) {
    markerGenes <- c(markerGenes, gene)
  }
}
print(markerGenes)

heatmapGI <- plotMarkerHeatmap(
  seMarker = markersGI,
  cutOff = "FDR <= 0.001 & Log2FC >= 0.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

markerList <- getMarkers(markersGI, cutOff = "FDR <= 0.01 & Log2FC >= 0.25")
for (cluster in selected_groups){
  write.csv(markerList[cluster], file= paste0(output_dir, "integration_marker_list_", cluster, ".csv"))
}
markerList$'B intermediate'
ComplexHeatmap::draw(heatmapGI,
                     heatmap_legend_side = "bot",
                     annotation_legend_side = "bot",
                     #row_order()
)
plotPDF(heatmapGI, name = "GeneIntegration-Marker-Heatmap",
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
  colorBy = "GeneIntegrationMatrix",
  name = markerGenes,
  embedding = "UMAP_n20_d0.1",
  imputeWeights = getImputeWeights(proj),
  continuousSet = "horizonExtra"
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

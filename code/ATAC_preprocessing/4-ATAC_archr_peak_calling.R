# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 10-31-2022
# Written by: Abhi Ramakrishnan
# Summary: Perform pseudo-bulking and peak calling for scATAC seq data
# Mem = 180G, Time = 100 hrs
# Important to run all pseudo-bulking and peakcalling (Step 1) all together
#-------------------------------------------------------------------------------
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
  library("BiocManager")
  library("rhdf5")
  library("BSgenome.Hsapiens.UCSC.hg38")
})

# Generate output directory
output_dir <- "path/to/output/directory/"

ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)

# Set random seed
set.seed(123)

# Source helper functions
source("path/to/helper_functions.R")

# Set number of threads
addArchRThreads(threads = 16)

# Load project
proj_dir <- "path/to/project/directory/"
proj <- loadArchRProject(proj_dir)

# Set working dir to project dir
setwd(proj_dir)

# Modify predictedGroup labels to remove spaces in cell type names
old_label <- proj$predictedGroup
new_label <- gsub(" ", "_", old_label)
table(old_label)
table(new_label)
proj$predictedGroupNew <- new_label
# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
#-------------------------------------------------------------------------------
# 1. Peak calling
#-------------------------------------------------------------------------------
# Find out number of cells per sample, per predicted RNA cell type
# Determine the parameters for pseudo-bulking based on this
print(table(proj$Sample, proj$predictedGroupNew))

# Perform pseudo-bulking of single cells
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "predictedGroupNew",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 500,
  maxFragments = 25 * 10^6,
  minReplicates = 3, # Changed from default 2
  maxReplicates = 52, # Changed from default 5 to total sample number
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE, # Changed this from false on 11.11.2022
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)

# Add peak set (NOTE: Each project can only contain one peak set!)
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "predictedGroupNew",
  pathToMacs2 = "macs2"

)
# Inspect peak set
getPeakSet(proj)

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
# Add peak matrix
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)


#-------------------------------------------------------------------------------
# 2. Identify marker peaks 
#-------------------------------------------------------------------------------

# Load project
proj <- loadArchRProject(proj_dir)

# Select only cell type groups with enough cells (remove ASDC and HSPC with 3 cells each)
predicted_groups <- data.frame(table(proj$predictedGroupNew))
selected_groups <- levels(predicted_groups$Var1)
selected_groups <- selected_groups[-c(1,19)]
# Identify marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroupNew",
  useGroups = selected_groups,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveHDF5SummarizedExperiment(markersPeaks, dir=output_dir, prefix="markersPeaks", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)
# Load markerPeaks object
markersPeaks <- loadHDF5SummarizedExperiment(dir=output_dir, prefix="markersPeaks")
markersPeaks
# Identify marker peaks for each cell type
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 0.125")
markerList
for (celltype in selected_groups) {
  write.csv(markerList[[celltype]], file= paste0(output_dir, "marker_list_", celltype, ".csv"))
}
# Generate marker heatpmap using peak matrix
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.125",
  pal = "blueYellow",
  clusterCols = TRUE,
  labelMarkers = NULL,
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

# Make volcano plots for differential peaks for each cell type cluster
celltypes_subset <- selected_groups[-c(8, 15, 17, 18, 24, 25)]
for (cellype in celltypes_subset) {
  pv <- plotMarkers(seMarker = markersPeaks,
                   name = celltype,
                   cutOff = "FDR <= 0.01 & abs(Log2FC) >= 0.5",
                   plotAs = "Volcano")
  set_panel_size(pv, dpi=300, file = paste0(output_dir, celltype, "_markers_volcano.pdf"),
                 width = unit(6, "in"), height = unit(6, "in"))
}


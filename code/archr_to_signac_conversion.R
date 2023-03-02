# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-25-23
# Written by: Natalie Piehl, Abhi Ramakrishnan
# Summary: Convert ArchR single cell ATAC object to Signac Seurat object,
# following the Swarup lab ArchRtoSignac tutorial
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load packages
library(ArchR)
library(Signac)
library(biovizBase)
library(Seurat)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(rlang)
library(devtools)
library(ArchRtoSignac)
library(rgeos)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)

# Generate output directory
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/conversion/archr_to_signac_cleaned_and_investigated/out_NP_01-26-2023/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)

# Set random seed
set.seed(123)

# Set threads
addArchRThreads(1)

# Define old proj dir
proj_dir <- "/projects/b1169/projects/AD_APOE/results/archr/project_main/"

# Define annotations path
annotations_path <- "/projects/b1169/projects/AD_APOE/results/archr/manual_outs/to_signac/2022_11_22_AR/UCSC.EnsDb.Hsapiens.v86.annotations.rds"

# Define new proj dirs
cd4_proj_dir <- "/projects/b1169/projects/AD_APOE/results/archr/project_cd4/"
noncd4_proj_dir <- "/projects/b1169/projects/AD_APOE/results/archr/project_noncd4/"

#-------------------------------------------------------------------------------
# Save cd4 and noncd4 ArchR projects
#-------------------------------------------------------------------------------

# Load ArchR project
proj <- loadArchRProject(proj_dir)

# Extract metadata from main proj
meta_sample_full <- data.frame(proj@sampleColData)
write.csv(meta_sample_full, paste0(output_dir, "full_archr_proj_sample_metadata.csv"))
meta_cell_full <- proj@cellColData
saveRDS(meta_cell_full, file = paste0(output_dir, "full_archr_proj_cell_metadata.rds"))

# Find CD4 T Cell barcodes
cd4_celltypes <- c("CD4_TCM", "CD4_TEM", "Treg",
                   "CD4_Naive", "CD4_Proliferating", "CD4_CTL")
cd4_ids <- BiocGenerics::which(proj$predictedGroupNew %in% cd4_celltypes)
cd4_cells <- proj$cellNames[cd4_ids]

# Subset for CD4 T Cells
subsetArchRProject(
  ArchRProj = proj,
  cells = cd4_cells,
  outputDirectory = cd4_proj_dir,
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)
proj

# Reload ArchR project
rm(proj)
proj <- loadArchRProject(proj_dir)

# Find non-CD4 T Cell barcodes
noncd4_ids <- BiocGenerics::which(proj$predictedGroupNew %!in% cd4_celltypes)
noncd4_cells <- proj$cellNames[noncd4_ids]

# Subset for non-CD4 T Cells
subsetArchRProject(
  ArchRProj = proj,
  cells = noncd4_cells,
  outputDirectory = noncd4_proj_dir,
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
proj

# Load in annotations
annotations <- readRDS(annotations_path)

#-------------------------------------------------------------------------------
# Convert ArchR to Signac  (CD4)
#-------------------------------------------------------------------------------

# the cellranger outputs directory
fragments_dir <- "/projects/b1042/Gate_Lab/AD_APOE/results/cellranger/atac/out_2022_08_11_NP/"

# Reload CD4 ArchR project
rm(proj)
proj <- loadArchRProject(cd4_proj_dir)

# Get peak matrix from ArchR object
cd4_pkm <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

# Generate Signac/Seurat object
cd4_s <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "hg38",
  fragments_dir = fragments_dir,
  fragments_fromcellranger = "Yes",
  pm = assay(cd4_pkm), # peak matrix from getPeakMatrix()
  annotation = annotations # annotation from getAnnotation()
)

# Extract UMAP
umap_df <- proj@embeddings$UMAP_n20_d0.1$df %>% as.matrix
rownames(umap_df) <- colnames(s) # make the rowname the same format as seurat
colnames(umap_df) <- c('UMAP_1', 'UMAP_2')

# Add to seurat object
cd4_s@reductions$umap <- Seurat::CreateDimReducObject(
  embeddings=umap_df,
  assay="peaks"
)

# Export object
saveRDS(cd4_s, file = paste0(output_dir, "cd4_s.rds"))

#-------------------------------------------------------------------------------
# Convert ArchR to Signac  (non-CD4)
#-------------------------------------------------------------------------------

# Reload nonCD4 ArchR project
rm(proj)
proj <- loadArchRProject(noncd4_proj_dir)

# Get peak matrix from ArchR object
noncd4_pkm <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

# Generate Signac/Seurat object
noncd4_s <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "hg38",
  fragments_dir = fragments_dir,
  fragments_fromcellranger = "Yes",
  pm = assay(noncd4_pkm), # peak matrix from getPeakMatrix()
  annotation = annotations # annotation from getAnnotation()
)

# Extract UMAP
umap_df <- proj@embeddings$UMAP_n20_d0.1$df %>% as.matrix
rownames(umap_df) <- colnames(s) # make the rowname the same format as seurat
colnames(umap_df) <- c('UMAP_1', 'UMAP_2')

# Add to seurat object
noncd4_s@reductions$umap <- Seurat::CreateDimReducObject(
  embeddings=umap_df,
  assay="peaks"
)

# Export object
saveRDS(noncd4_s, file = paste0(output_dir, "noncd4_s.rds"))
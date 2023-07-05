# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 02-06-2023
# Written by: Natalie Piehl
# Summary: Run TF-IDF normalization
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
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
})

# Organize inputs
noncd4_s_path <- "/path/to/noncd4_seurat_object"
cd4_s_path <- "/path/to/cd4_seurat_object"
output_dir <- "/path/to/output_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run TF-IDF

# Load Seurat object
s <- readRDS(noncd4_s_path)
s

# Run TF-IDF
s <- RunTFIDF(s)
s

# Export result
saveRDS(s, file = paste0(output_dir, "noncd4_s_TFIDF.rds"))

# Load Seurat object
s <- readRDS(cd4_s_path)
s

# Run TF-IDF
s <- RunTFIDF(s)
s

# Export result
saveRDS(s, file = paste0(output_dir, "cd4_s_TFIDF.rds"))

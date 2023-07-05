# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-06-2022
# Written by: Natalie Piehl
# Summary: Prep pseudobulk coverages from ArchR
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("ArchR")
  library("data.table")
})

# Organize inputs
proj_dir <- "/path/to/ArchR/proj/"
cellranger_dir <- "/path/to/cellranger/results/"
output_dir <- "/path/to/output_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set threads
addArchRThreads(1)

#-------------------------------------------------------------------------------
# Export fragments
#-------------------------------------------------------------------------------

# Load ArchR project
proj <- loadArchRProject(proj_dir)

# Add broad celltype column
proj@cellColData$broad_celltype <- mapvalues(as.vector(proj@cellColData$predictedGroupNew),
                                        from = c("B_intermediate", "B_memory", "B_naive", "Plasmablast",
                                                 "CD14_Mono", "CD16_Mono",
                                                 "ASDC", "cDC1", "cDC2", "pDC",
                                                 "CD4_CTL", "CD4_Naive", "CD4_Proliferating", "CD4_TCM", "CD4_TEM", "Treg",
                                                 "CD8_Naive", "CD8_Proliferating", "CD8_TCM", "CD8_TEM", "MAIT",
                                                 "NK", "NK_Proliferating", "NK_CD56bright",
                                                 "ILC", "dnT", "gdT",
                                                 "Platelet", "Eryth", "HSPC", "Doublet"),
                                        to = c(rep("B_Cells", 4),
                                               rep("Monocytes", 2),
                                               rep("Dendritic_Cells", 4),
                                               rep("CD4+_T_Cells", 6),
                                               rep("CD8+_T_Cells", 5),
                                               rep("NK_Cells", 3),
                                               rep("Other_T_Cells", 3),
                                               rep("Other", 4)))

# Make sample+celltype column for grouping
proj@cellColData$sample_celltype <- paste(proj@cellColData$Sample, proj@cellColData$broad_celltype, sep = "_")

# Iterate through each sample+celltype column
for (sample in unique(proj@cellColData$Sample)) {
  # Read in fragments
  frags <- fread(paste0(cellranger_dir, sample, "/outs/fragments.tsv.gz"))
  names(frags) <- c("chrom", "chromStart", "chromEnd", "barcode", "readSupport")
  
  for (celltype in unique(proj@cellColData$broad_celltype)) {
    # Define group
    group <- paste(sample, celltype, sep = "_")
  
    # Extract barcodes
    barcodes <- rownames(proj@cellColData[which(proj@cellColData$sample_celltype == group),])
    barcodes <- sapply(barcodes, function(x) {gsub(".*#", "", x)}) %>% as.vector
    
    # Extract fragments of interest
    frags_oi <- frags[which(frags$barcode %in% barcodes),]
  
    # Export as tsv.gz
    write.table(frags_oi, paste0(output_dir, group, "_fragments.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}
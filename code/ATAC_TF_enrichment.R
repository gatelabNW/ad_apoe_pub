# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-10-2023
# Written by: Natalie Piehl
# Summary: Run TF enrichment on DARs
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
  library("JASPAR2020")
  library("TFBSTools")
  library("BSgenome.Hsapiens.UCSC.hg38")
})

# Organize inputs
celltype_colors_path <- "/path/to/celltype_color_map.csv"
da_base_dir <- "/path/to/da/results/"
output_base_dir <- "/path/to/output/folder/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# List comparisons
comparisons <- c("ADvsHC",
                 "ADvsHC_33", "ADvsHC_34", "ADvsHC_44",
                 "34vs33_HC", "34vs33_AD",
                 "44vs33_HC", "44vs33_AD",
                 "44vs34_HC", "44vs34_AD")

# List types
types <- c("noncd4", "cd4")

#-------------------------------------------------------------------------------
# Run TF enrichment
#-------------------------------------------------------------------------------

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Define celltypes
color_map <- read.csv(celltype_colors_path)
cell_types <- sapply(color_map$predicted.celltype.l2,
                     function(x) {gsub(" ", "_", x)}) %>% as.vector

for (type in types) {
  # Load Seurat object
  s <- readRDS("path/to/type/seurat/object/")
  
  # add motif information
  s <- AddMotifs(s, 
                 genome = BSgenome.Hsapiens.UCSC.hg38,
                 pfm = pwm)
  
  # Add broad cell types
  s@meta.data$broad_celltype <- mapvalues(s@meta.data$predictedGroupNew,
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
  
  for (comparison in comparisons[2:4]) {
    # Define comparison output dir
    output_dir <- paste0(output_base_dir, comparison, "/")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    run_enrichment <- function(cell_type) {
      if (!file.exists(paste0(output_dir, cell_type, "_", comparison, "_enriched_motifs_inUpDARs.csv"))) {
        tryCatch({
          # Read in DARs
          print(cell_type)
          da_peaks_path <- paste0(da_base_dir, comparison, "/", cell_type, "_", comparison, "_dars.csv")
          da_peaks <- read.csv(da_peaks_path)
          
          # Subset seurat object
          tmp <- subset(s, broad_celltype == cell_type)
          
          #---------------------------------------------------------------------
          # Enriched in Up DARs
          print("Analyzing TFs enriched in Up DARs...")
          
          # Identify up-regulated DARs
          top_da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & da_peaks$avg_log2FC > lfc.thresh), "sitename"] %>% as.vector
          
          if (length(top_da_peaks) == 0) {
            print(paste0("No up DARs in ", cell_type))
          } else {
            print(paste(as.character(length(top_da_peaks)), " up DARs"))
            
            # Test enrichment
            enriched_motifs <- FindMotifs(
              object = tmp,
              features = top_da_peaks,
              assay = "peaks"
            )
            
            # Export enriched motifs
            write.csv(enriched_motifs, file = paste0(output_dir, cell_type, "_", comparison, "_enriched_motifs_inUpDARs.csv"))
          }
          
          #---------------------------------------------------------------------
          # Enriched in Down DARs
          print("Analyzing TFs enriched in Down DARs...")
          
          # Identify down-regulated DARs
          top_da_peaks <- da_peaks[which(da_peaks$BH < padj.thresh & da_peaks$avg_log2FC < -lfc.thresh), "sitename"] %>% as.vector
          
          if (length(top_da_peaks) == 0) {
            print(paste0("No down DARs in ", cell_type))
          } else {
            print(paste(as.character(length(top_da_peaks)), " down DARs"))
            
            # Test enrichment
            enriched_motifs <- FindMotifs(
              object = tmp,
              features = top_da_peaks,
              assay = "peaks"
            )
            
            # Export enriched motifs
            write.csv(enriched_motifs, file = paste0(output_dir, cell_type, "_", comparison, "_enriched_motifs_inDownDARs.csv"))
          }
          
          # Remove seurat object
          rm(tmp)
        }, error = function(e) {
          message(e)
        })  
      }
    }
    lapply(cell_types, run_enrichment)
  }
}
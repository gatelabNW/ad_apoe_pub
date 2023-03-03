# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-11-2022
# Written by: Natalie Piehl
# Summary: Remove background gene expression and export corrected counts via SoupX
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Matrix")
  library("SoupX")
  library("DropletUtils")
  library("purrr")
  library("Hmisc")
  library("doMC")
})

# Initialize input parameters
cellranger_dir <- "/projects/b1042/Gate_Lab/AD_APOE/results/cellranger/rna/out_2022_08_10_NP"
output_dir <- "/projects/b1169/projects/AD_APOE/results/soupx/main/out_NP_08-16-2022/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
plot_dir <- paste0(output_dir, "qc/")
dir.create(plot_dir, showWarnings = FALSE)

# Set core number for parallel model fitting
registerDoMC(cores = 12)

# Create lists of directories to load as SoupX objects and Seurat objects
sample_dirs <- list.dirs(cellranger_dir, recursive = FALSE)
print(sample_dirs)

# Initialize gene list to estimate contamination fraction
# Selecting for monocyte/dendritic markers which are highly specific and abundant
M.genes <- list(MC = c("S100A8", "S100A9"))

# For each sample...
run_soupx <- function(dir) {
  # Isolate sample name
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)

  # Make sample output dir
  sample_dir <- paste0(output_dir, "/", sample, "/")
  dir.create(sample_dir, showWarnings = FALSE)

  # Print what sample is being processed
  print(paste0("Processing sample ", sample))

  # Load in SoupX and Seurat data
  soupx <- load10X(paste0(dir, "/outs"))
  seurat <- Read10X(paste0(dir, "/outs/filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  seurat

  # Normalize and run PCA on Seurat object
  seurat <- SCTransform(seurat, variable.features.n = 2000, verbose = TRUE) %>%
    RunPCA()

  # Generate clusters
  seurat <- RunTSNE(seurat, dims = 1:10) %>%
    FindNeighbors(dims = 1:10) %>%
    FindClusters(resolution = 0.3)

  # Add cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))

  # Set fraction manually for failed samples
  manual_samples <- c("G1010_y2", "G1028", "G1034", "G738", "G773_y4")
  if (sample %in% manual_samples) {
    soupx <- setContaminationFraction(soupx, 0.01425)
  } else {
    # Evaluate marker genes
    write.csv(head(soupx$soupProfile[order(soupx$soupProfile$est, decreasing = TRUE), ], n = 50),
              file = paste0(plot_dir, sample, "_top_background_genes.csv"))
    p <- plotMarkerDistribution(soupx)
    set_panel_size(p, file = paste0(plot_dir, sample, "_marker_plot.pdf"))

    # Estimate contamination fraction
    useToEst <- estimateNonExpressingCells(soupx, nonExpressedGeneList = M.genes)
    soupx <- calculateContaminationFraction(soupx,
                                            M.genes,
                                            useToEst = useToEst,
                                            forceAccept = TRUE)
  }

  # Create adjusted counts
  print(paste0("Fraction: ", sample, " ", soupx$metaData$rho[1]))
  adj_counts <- adjustCounts(soupx)

  # Export change map
  p <- plotChangeMap(soupx, adj_counts, "S100A8")
  set_panel_size(p, file = paste0(plot_dir, sample, "_S100A8_changeplot.pdf"))
  p <- plotChangeMap(soupx, adj_counts, "S100A9")
  set_panel_size(p, file = paste0(plot_dir, sample, "_S100A9_changeplot.pdf"))

  # Save Corrected Counts
  DropletUtils:::write10xCounts(paste0(output_dir, "/", sample), adj_counts, overwrite = TRUE)

  return(soupx$metaData$rho[1])
}

# Run SoupX
rho <- mclapply(sample_dirs, run_soupx, mc.cores = 12)
rho

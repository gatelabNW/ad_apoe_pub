# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-16-2022
# Written by: Natalie Piehl
# Summary: Preprocess RNA with Seurat
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
})

# Initialize input parameters
soupx_dir <- "/projects/b1169/projects/AD_APOE/results/soupx/main/out_NP_08-16-2022/"
metadata_path <- "/projects/b1169/projects/AD_APOE/data/metadata/final_samples_AD_APOE.xlsx"
output_dir <- "/projects/b1169/projects/AD_APOE/results/seurat/rna_preprocessing/out_NP_08-19-2022_no1086/"
input_dir <- "/projects/b1169/projects/AD_APOE/results/seurat/rna_preprocessing/out_NP_08-16-2022/"
prefiltering_qc_plot_dir <- paste0(output_dir, "prefiltering_qc_plots/")
postfiltering_qc_plot_dir <- paste0(output_dir, "postfiltering_qc_plots/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(prefiltering_qc_plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(postfiltering_qc_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Set core number for parallel model fitting
registerDoMC(cores = 12)

#-------------------------------------------------------------------------------
# Generate seurat objects

# List SoupX corrected counts sample dirs
soupx_sample_dirs <- list.dirs(path = soupx_dir, full.names = TRUE, recursive = FALSE)
samples <- list.dirs(path = soupx_dir, full.names = FALSE, recursive = FALSE)
samples

# Create function to load seurat objects
load_seurat <- function(dir) {
  counts <- Read10X(dir)
  project <- tail(unlist(strsplit(dir, "/")), 1)
  return(CreateSeuratObject(counts = counts, project = project,
                            min.cells = 3, min.features = 200))
}

# Apply Seurat loading function to all samples
seurat_object_list <- sapply(soupx_sample_dirs, load_seurat)

#-------------------------------------------------------------------------------
# Identify doublets with DoubletFinder

# Define doublet formation rate
doublet_formation_rate <- 0.054
print(paste0("Using doublet formation rate of ", doublet_formation_rate))

run_doubletfinder <- function(s) {
  ## Process normally
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Run TSNE clustering
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunTSNE(s, dims = 1:12)

  ## pK Identification (no ground-truth) ---------------------------------------
  sweep.res.list <- paramSweep_v3(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  ## Homotypic Doublet Proportion Estimate -------------------------------------
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print(s)

  ## Run DoubletFinder with varying classification stringencies ----------------
  s <- doubletFinder_v3(s, PCs = 1:12, pN = 0.25, pK = optimal_pK,
                        nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  ## Rename column name for consistency
  colnames(s@meta.data)[ grep("DF.classifications*",
                              colnames(s@meta.data)) ] <- "DF.classifications"
  print(table(s[["DF.classifications"]]))

  return(s)
}

# Apply DoubletFinder function to all samples
seurat_object_list <- mclapply(seurat_object_list, run_doubletfinder,
                               mc.cores = 12)

#-------------------------------------------------------------------------------
# Merge into one seurat object

# Merge objects
s <- merge(seurat_object_list[[1]],
           unlist(seurat_object_list,
                  use.names = FALSE)[2:length(seurat_object_list)],
           add.cell.ids = samples,
           project = "AD_APOE")
s

# Save result
save(s, file = paste0(output_dir, "s_rawmerge"))
# load(paste0(input_dir, "s_rawmerge"))

#-------------------------------------------------------------------------------
# Add Metadata

# Load in metadata
meta <- read.xlsx(metadata_path, sheetIndex = 3)

# Extract cell barcodes and IDs
ids <- s[["orig.ident"]]

# Add formatted id to metadata
meta$id <- paste0("G", meta$Patient_ID)
for (i in seq(nrow(meta))) {
  if (meta$Study_year[i] != 1) {
    meta$id[i] <- paste0(meta$id[i], "_y", meta$Study_year[i])
  }
}
meta[which(meta$id %in% c("G1241_y3", "G906_y3")),
     "id"] <- c("G1241_Y3", "G906_Y3")
meta <- relocate(meta, id)

# Subset for metadata columns of interest
meta <- meta[,c("id", "Age", "Sex", "Diagnosis",
                "APOE_genotype", "Used_for_exp")]
names(meta) <- c("id", "Age", "Sex", "Diagnosis",
                 "APOE_genotype", "Experiment")

# Map experiments
meta$Experiment <- mapvalues(meta$Experiment,
                             from = c("", "N", "Y"),
                             to = c("RNA_ATAC_diffyear",
                                    "RNA",
                                    "RNA_ATAC"))

# Add on sort day
meta$Sort_day <- c(rep(1,3), rep(2,4), rep(3,3), rep(4,7),
                   rep(5,7), rep(6,7), rep(7, 8), rep(8, 6), 3, rep(9,6))

# Remove all pANN columns
s@meta.data <- s@meta.data %>% select(-contains("pANN"))
s[["RNA_snn_res.0.3"]] <- NULL
s[["seurat_clusters"]] <- NULL

# Merge metadata onto cell barcodes
meta_bycell <- merge(ids, meta, all.x = TRUE,
                     by.x = "orig.ident", by.y = "id")
rownames(meta_bycell) <- rownames(ids)

# Add on metadata
s <- AddMetaData(object = s,
                 metadata = meta_bycell)

# Set levels
s@meta.data$Diagnosis <- factor(s@meta.data$Diagnosis,
                                levels = c("Healthy Control", "Alzheimers Disease"))
s@meta.data$Age <- factor(s@meta.data$Age,
                          levels = seq(100))
s@meta.data$Sort_day <- factor(s@meta.data$Sort_day,
                               levels = seq(10))

#-------------------------------------------------------------------------------
# Generate QC plots

# Calculate Percent Mitochondrial Genes
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")

run_qc <- function(qc_plot_dir) {
  # Cell count per sample
  Idents(s) <- 'orig.ident'
  cell_counts <- table(s[["orig.ident"]]) %>% as.data.frame
  myplot <- ggplot(cell_counts, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat="identity") +
    theme_Publication_blank() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sample ID",
         y = "Number of Doublets",
         title = "Number of Cells per Sample ID")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCells_ID.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  # Doublet count per sample
  doublet_counts <- table(s[[c("orig.ident", "DF.classifications")]]) %>% as.data.frame
  myplot <- ggplot(doublet_counts, aes(x = orig.ident, y = Freq, fill = DF.classifications)) +
    geom_bar(stat="identity") +
    theme_Publication_blank() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sample ID",
         y = "Number of Doublets",
         title = "Number of Doublets per Sample ID")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nDoublets_ID.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  # By ID
  Idents(s) <- 'orig.ident'
  myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sample ID",
         y = "Number of Features",
         title = "Number of Features per Sample ID")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_ID.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sample ID",
         y = "Number of Counts",
         title = "Number of Counts per Sample ID")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_ID.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sample ID",
         y = "Percent Mitochondrial",
         title = "Percent Mitochondrial per Sample ID")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_ID.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  # By Age
  Idents(s) <- 'Age'
  myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Number of Features (Genes)",
         title = "nFeatures per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Age.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Number of Counts",
         title = "nCounts per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Age.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Percent Mitochondrial",
         title = "Percent MT per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Age.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  # By Diagnosis
  Idents(s) <- 'Diagnosis'
  myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0, cols = c("black", "red")) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Diagnosis",
         y = "Number of Features",
         title = "Number of Features per Diganosis")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_Diagnosis.pdf"),
                 width=unit(3, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0, cols = c("black", "red")) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Diagnosis",
         y = "Counts",
         title = "Number of Counts per Diganosis")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_Diagnosis.pdf"),
                 width=unit(3, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0, cols = c("black", "red")) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Diagnosis",
         y = "Percent Mitochondrial",
         title = "Percent Mitochondrial per Diagnosis")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_Diagnosis.pdf"),
                 width=unit(3, "in"), height=unit(4, "in"))

  # By sort day
  Idents(s) <- 'Sort_day'
  myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sort Day",
         y = "Number of Features (Genes)",
         title = "nFeatures per Sort Day")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_sort_day.pdf"),
                 width=unit(4, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sort Day",
         y = "Number of Counts",
         title = "nCounts per Sort Day")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_sort_day.pdf"),
                 width=unit(4, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Sort Day",
         y = "Percent Mitochondrial",
         title = "Percent MT per Sort Day")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_sort_day.pdf"),
                 width=unit(4, "in"), height=unit(4, "in"))

  # By Doublets
  Idents(s) <- 'DF.classifications'
  myplot <- VlnPlot(s, features = c("nFeature_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Number of Features (Genes)",
         title = "nFeatures per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_doublets.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("nCount_RNA"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Number of Counts",
         title = "nCounts per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nCount_doublets.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  myplot <- VlnPlot(s, features = c("percent.mt"), pt.size = 0) +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Age",
         y = "Percent Mitochondrial",
         title = "Percent MT per Age")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "percent.mt_doublets.pdf"),
                 width=unit(8, "in"), height=unit(4, "in"))

  # Look at Counts vs Features
  Idents(s) <- 'orig.ident'
  myplot <- FeatureScatter(s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    theme_Publication_blank() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 20000, linetype = 2, color = "gray") +
    geom_hline(yintercept = 4000, linetype = 2, color = "gray") +
    labs(x = "Number of Counts",
         y = "Number of Features (Genes)",
         title = "nFeatures per nCounts")
  set_panel_size(myplot, file=paste0(qc_plot_dir, "nFeature_by_nCount.pdf"),
                 width=unit(6, "in"), height=unit(4, "in"))
}

# Run QC
run_qc(prefiltering_qc_plot_dir)

#-------------------------------------------------------------------------------
# Filtering and rerun QC

# Subset for singlets
s <- subset(s, subset = DF.classifications == "Singlet")

print("Before removing high mitochondrial content cells...")
s
# Remove cells with above threshold mitochondrial reads
s <- subset(s, subset = percent.mt < 10)
print("After removing high mitochondrial content cells...")
s

# Run QC postfilter
run_qc(postfiltering_qc_plot_dir)

#-------------------------------------------------------------------------------
# Run naive dimension reduction

# Generate PCA with batch effects
s <- NormalizeData(s,
                   normalization.method = "LogNormalize",
                   scale.factor = 10000)
s <- FindVariableFeatures(s,
                          selection.method = "vst",
                          nfeatures = 2000)
s <- ScaleData(s)
s <- RunPCA(object = s)

# Cluster Cells
s <- FindNeighbors(s, dims = 1:12)
s <- FindClusters(s, resolution = 0.3)
s <- RunUMAP(s, dims = 1:12)

# Make binned counts, features, and mt percent
s@meta.data$count_bin <- cut(s@meta.data$nCount_RNA, 10)
s@meta.data$feature_bin <- cut(s@meta.data$nFeature_RNA, 10)
s@meta.data$percentmt_bin <- cut(s@meta.data$percent.mt, 10)

# PCA labelled by sort day
s <- SetIdent(s, value = 'Sort_day')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_sort_day.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by ID
s <- SetIdent(s, value = 'orig.ident')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_ID.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by diagnosis
s <- SetIdent(s, value = 'Diagnosis')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0, cols = c("black", "red"))
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_diagnosis.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by counts
s <- SetIdent(s, value = 'count_bin')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_nCount.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by features
s <- SetIdent(s, value = 'feature_bin')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_nFeature.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by mitochondrial
s <- SetIdent(s, value = 'percentmt_bin')
myplot <- DimPlot(s, reduction = "pca", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "PCA_percentmt.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# UMAP labelled by sort day
s <- SetIdent(s, value = 'Sort_day')
myplot <- DimPlot(s, reduction = "umap", pt.size = 1, shuffle = TRUE)
set_panel_size(myplot, file = paste0(postfiltering_qc_plot_dir, "UMAP_sort_day.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# UMAP labelled by ID
s <- SetIdent(s, value = 'orig.ident')
myplot <- DimPlot(s, reduction = "umap", pt.size = 0)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "UMAP_ID.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# UMAP labelled by diagnosis
s <- SetIdent(s, value = 'Diagnosis')
myplot <- DimPlot(s, reduction = "umap", pt.size = 0, cols = c("black", "red"))
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "UMAP_diagnosis.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# UMAP labelled by counts
s <- SetIdent(s, value = 'count_bin')
myplot <- DimPlot(s, reduction = "umap", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "UMAP_nCount.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# UMAP labelled by features
s <- SetIdent(s, value = 'feature_bin')
myplot <- DimPlot(s, reduction = "umap", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "UMAP_nFeature.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# PCA labelled by mitochondrial
s <- SetIdent(s, value = 'percentmt_bin')
myplot <- DimPlot(s, reduction = "umap", pt.size = 0, shuffle = TRUE)
set_panel_size(myplot, file=paste0(postfiltering_qc_plot_dir, "UMAP_percentmt.pdf"),
               width=unit(5, "in"), height=unit(4, "in"))

# Export seurat object
save(s, file = paste0(output_dir, "s_qced"))
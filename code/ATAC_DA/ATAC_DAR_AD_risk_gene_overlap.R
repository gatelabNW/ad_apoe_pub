# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-04-2023
# Written by: Natalie Piehl
# Summary: Investigate AD risk genes in DARs
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("data.table")
  library("pheatmap")
  library("corrplot")
  library("GenomicRanges")
})

# Organize inputs
da_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/"
intersection_dir <- "/projects/b1169/projects/AD_APOE/results_atac/da_broad_celltypes/LR_DESeq2/out_NP_05-22-2023/full_intersection/"
ranges_path <- "/projects/b1169/projects/AD_APOE/data/ranges/full_ranges.rds"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/ad_risk_genes/dars/out_NP_05-22-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set thresholds
lfc.thresh <- .125
padj.thresh <- .05

# Specify comparison
comparisons <- c("ADvsHC", "ADvsHC_44")

#-------------------------------------------------------------------------------
# Export AD risk DARs

# Define risk genes in order
risk_genes <- c("APP", "PSEN1", "PSEN2",
                "APOE", "BIN1", "MS4A6A", "PICALM", "CR1", "CLU", "TREM2",
                "ABCA7", "NYAP1", "PTK2B", "PLCG2", "SPI1", "SORL1", "HLA-DRB1",
                "CD2AP", "SLC24A4", "RIN3", "ADAMTS1", "ADAMTS4", "CASS4",
                "ADAM10", "FERMT2", "HAVCR2", "SCIMP", "CLNK", "ECHDC3",
                "TNIP1", "ABCA1", "CNTNAP2", "USP6NL", "INPP5D", "CD33", "ACE",
                "IQCK", "WWOX", "ABI3", "HESX1", "FHL2", "APH1B", "HS3ST1",
                "CHRNE", "CCDC6", "AGRN", "KAT8", "IL34")

# Define cell types
cell_types <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells", "Other_T_Cells")

for (comparison in comparisons) {
  # Read in shared DARs
  shared_dars <- read.csv(paste0(intersection_dir, comparison, "_padj0.05_lfc0.125_LR+DESeq2_DARs.csv"))
  
  # Define comparison deg dirs
  if (comparison %in% c("ADvsHC_44", "44vs34_AD")) {
    dar_dir <- paste0(da_base_dir, "main/out_NP_02-06-2023/", comparison, "/")
  } else {
    dar_dir <- paste0("/projects/b1169/projects/AD_APOE/results_atac/LR_withEthnicity/out_NP_05-22-2023/", comparison, "/")
  }

  # Iterate through cell types
  find_risk_gene_dars <- function(cell_type) {
    tryCatch({
      if (cell_type == "CD8+_T_Cells") {
        celltype_files <- list.files(dar_dir, pattern = "CD8", full.names = TRUE)
      } else if (cell_type == "CD4+_T_Cells") {
        celltype_files <- list.files(dar_dir, pattern = "CD4", full.names = TRUE)
      } else {
        celltype_files <- list.files(dar_dir, pattern = cell_type, full.names = TRUE)
      }
      da_peaks <- read.csv(grep(".csv", celltype_files, value = TRUE))
      
      # Subset for shared with DESeq2
      da_peaks <- da_peaks[which(da_peaks$sitename %in% shared_dars[,gsub("\\+", ".", cell_type)]),]
      
      # Subset for sig peaks
      sig_peaks <- da_peaks[which(da_peaks$BH <= padj.thresh & abs(da_peaks$avg_log2FC) >= lfc.thresh),]
      
      # Isolate peaks associate with risk genes
      sig_peaks <- sig_peaks[which(sig_peaks$nearestGene %in% risk_genes),]
      sig_peaks$cell_type <- rep(cell_type, nrow(sig_peaks))
      
      return(sig_peaks)
    }, error = function(e) {
      message(paste(cell_type, e))
    })
  }
  
  # Run iteration
  res <- lapply(cell_types, find_risk_gene_dars)
  
  # Combine
  res_df <- rbindlist(res)
  
  # Export result
  write.csv(res_df, paste0(output_dir, "ad_risk_genes_in_", comparison, "_DARs.csv"))
  
  # Iterate through cell types
  num_risk_gene_dars <- function(cell_type) {
    tryCatch({
      if (cell_type == "CD8+_T_Cells") {
        celltype_files <- list.files(dar_dir, pattern = "CD8", full.names = TRUE)
      } else if (cell_type == "CD4+_T_Cells") {
        celltype_files <- list.files(dar_dir, pattern = "CD4", full.names = TRUE)
      } else {
        celltype_files <- list.files(dar_dir, pattern = cell_type, full.names = TRUE)
      }
      da_peaks <- read.csv(grep(".csv", celltype_files, value = TRUE))
      
      # Subset for shared with DESeq2
      da_peaks <- da_peaks[which(da_peaks$sitename %in% shared_dars[,gsub("\\+", ".", cell_type)]),]
      
      # Subset for sig peaks
      sig_peaks <- da_peaks[which(da_peaks$BH <= padj.thresh & abs(da_peaks$avg_log2FC) >= lfc.thresh),]
      
      # Isolate peaks associate with risk genes
      sig_peaks <- sig_peaks[which(sig_peaks$nearestGene %in% risk_genes),]
      
      # Generate gene dar frequencies
      risk_gene_freq <- table(sig_peaks$nearestGene)
      
      # Create formatted
      risk_gene_freq_formatted <- list()
      for (risk_gene in risk_genes) {
        if (risk_gene %in% names(risk_gene_freq)) {
          risk_gene_freq_formatted[[risk_gene]] <- risk_gene_freq[[risk_gene]]
        } else {
          risk_gene_freq_formatted[[risk_gene]] <- 0
        }
      }
      
      # Add cell type
      risk_gene_freq_formatted[["cell_type"]] <- cell_type
      
      return(risk_gene_freq_formatted)
    }, error = function(e) {
      message(paste(cell_type, e))
    })
  }
  
  # Run iteration
  res <- lapply(cell_types, num_risk_gene_dars)
  
  # Combine
  res_df <- rbindlist(res) %>% as.data.frame
  
  # Convert rownames to cell_type
  rownames(res_df) <- res_df$cell_type
  res_df$cell_type <- NULL
  
  # Generate with hclust
  data_heatmap <- t(res_df)
  colfunc <- colorRampPalette(c("white", "purple"))
  pheatmap(as.matrix(data_heatmap),
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           color = colfunc(5),
           breaks=seq(0, 4, length.out = 6),
           # border_color = NA,
           filename = paste0(output_dir, "num_", comparison, "_DARs_in_ad_risk_genes_heatmap.pdf"),
           width = 3,
           height = 7
  )
}

#-------------------------------------------------------------------------------
# Visualize AD risk gene and celltype at a time

# Define gene and comparison
gene <- "BIN1"
comparison <- "ADvsHC_44"

# Read in ranges and isolate gene OI
ranges <- readRDS(ranges_path)
ranges <- ranges[which(ranges$nearestGene == gene),]

# Load Seurat objects
cd4_s <- readRDS("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds")
cd4_subset <- subset(cd4_s, features = ranges$sitename)
rm(cd4_s)
noncd4_s <- readRDS("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds")
noncd4_subset <- subset(noncd4_s, features = ranges$sitename)
rm(noncd4_s)
s <- merge(cd4_subset, noncd4_subset)
rm(noncd4_subset)
rm(cd4_subset)

# Add tmp column
s@meta.data$tmp <- rep("tmp", nrow(s[[]]))

# Find min and max bps
gene_ranges <- granges(s) %>% data.frame
min_bp <- min(gene_ranges$start)
max_bp <- max(gene_ranges$end)

# Define gene range
range_oi <- GRanges(seqnames=unique(gene_ranges$seqnames), ranges=paste(min_bp, max_bp, sep = "-"))

# Set identity
Idents(s) <- "tmp"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = range_oi,
  extend.downstream = 50,
  extend.upstream = 50,
  annotation = TRUE,
  peaks = TRUE
)

# Export plot
pdf(paste0(output_dir, gene, "_tracks.pdf"),
    width = 8, height = 3)
print(cov_plot)
dev.off()

# # Write out gene ranges
# write.csv(gene_ranges, file = paste0(output_dir, gene, "_ranges.csv"))

# Read in shared DARs
shared_dars <- read.csv(paste0(intersection_dir, comparison, "_padj0.05_lfc0.125_LR+DESeq2_DARs.csv"), row.names = 1)

# Define comparison deg dirs
if (comparison %in% c("ADvsHC_44", "44vs34_AD")) {
  dar_dir <- paste0(da_base_dir, "main/out_NP_02-06-2023/", comparison, "/")
} else {
  dar_dir <- paste0("/projects/b1169/projects/AD_APOE/results_atac/LR_withEthnicity/out_NP_05-22-2023/", comparison, "/")
}

# Iterate through cell types
specific_risk_gene_dars <- function(cell_type) {
  tryCatch({
    if (cell_type == "CD8._T_Cells") {
      celltype_files <- list.files(dar_dir, pattern = "CD8", full.names = TRUE)
    } else if (cell_type == "CD4._T_Cells") {
      celltype_files <- list.files(dar_dir, pattern = "CD4", full.names = TRUE)
    } else {
      celltype_files <- list.files(dar_dir, pattern = cell_type, full.names = TRUE)
    }
    da_peaks <- read.csv(grep(".csv", celltype_files, value = TRUE))
    
    # Subset for shared with DESeq2
    da_peaks <- da_peaks[which(da_peaks$sitename %in% shared_dars[,cell_type]),]
    
    # Isolate sig peaks
    sig_peaks <- shared_dars[,gsub("\\+", ".", cell_type)]
    
    # Annotate with ranges
    sig_peaks <- da_peaks[which(da_peaks$sitename %in% ranges$sitename),]
    
    # Add cell type
    sig_peaks$cell_type <- rep(cell_type, nrow(sig_peaks))
    
    return(sig_peaks)
  }, error = function(e) {
    message(e)
  })
}

# Run func
res <- lapply(names(shared_dars), specific_risk_gene_dars)

# Combine
res_df <- rbindlist(res) %>% as.data.frame
res_df <- res_df[order(res_df$sitename),]

# Generate coverage plot highlighting DARs
gene_dars <- unique(res_df$sitename)
coords <- sapply(gene_dars, function(x) {gsub(paste0(unique(gene_ranges$seqnames), "-"), "", x)}) %>% unlist() %>% as.vector()
highlight <- GRanges(seqnames=unique(gene_ranges$seqnames), ranges=coords)
cov_plot <- CoveragePlot(
  object = s,
  region = range_oi,
  region.highlight = highlight,
  extend.downstream = 50,
  extend.upstream = 50,
  annotation = TRUE,
  peaks = TRUE
)

# Export plot
pdf(paste0(output_dir, gene, "_tracks_", comparison, "_dar_highlight.pdf"),
    width = 8, height = 3)
print(cov_plot)
dev.off()

# Isolate lfc values and padj values
r_mat <- res_df[,c("sitename", "cell_type", "avg_log2FC")]
p_mat <- res_df[,c("sitename", "cell_type", "BH")]

# Convert to wide (cell type row, peak col)
r_mat_wide <- pivot_wider(r_mat,
                          id_cols = cell_type,
                          values_from = avg_log2FC,
                          names_from = sitename)
p_mat_wide <- pivot_wider(p_mat,
                          id_cols = cell_type,
                          values_from = BH,
                          names_from = sitename)

# Convert NAs
r_mat_wide[is.na(r_mat_wide)] <- 0
p_mat_wide[is.na(p_mat_wide)] <- 1

# # for BIN1 with over 1 LFC, adjust legend in post
# r_mat_wide[,2:ncol(r_mat_wide)] <- r_mat_wide[,2:ncol(r_mat_wide)] / 1.25

# Convert to data.frame
r_mat_wide <- data.frame(r_mat_wide)
p_mat_wide <- data.frame(p_mat_wide)

# Remove cell_type
rownames(r_mat_wide) <- r_mat_wide$cell_type
r_mat_wide$cell_type <- NULL
r_mat_wide <- as.matrix(r_mat_wide)
rownames(p_mat_wide) <- p_mat_wide$cell_type
p_mat_wide$cell_type <- NULL
p_mat_wide <- as.matrix(p_mat_wide)

# Generate plot
pdf(paste0(output_dir, gene, "_in_", comparison, "_DARs_heatmap.pdf"),
    width = 4, height = 8)
corrplot(r_mat_wide,
         tl.col = "black",
         p.mat = p_mat_wide,
         method = "color",
         col=colorRampPalette(c("blue","white","red"))(200))
text(col(p_mat_wide), rev(row(p_mat_wide)),
     formatC(p_mat_wide, format = "e", digits = 2),
     cex=0.5)
dev.off()
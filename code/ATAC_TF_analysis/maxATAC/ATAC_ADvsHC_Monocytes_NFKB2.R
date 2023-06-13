# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-18-2022
# Written by: Natalie Piehl
# Summary: investgigate maxatac results
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
  library("ggbreak")
  library("ChIPseeker")
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  library("rtracklayer")
  library("org.Hs.eg.db")
})

# Organize inputs
res_dir <- "/projects/b1169/projects/AD_APOE/results_atac/maxatac/peaks/out_NP_04-25-2023/"
noncd4_seurat_object <- "/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/maxatac/NFKB2/out_NP_06-07-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparisons of interest
comparisons_oi <- c("HC_Monocytes", "AD_Monocytes")

# Define cell type of interest
cell_type <- "Monocytes"

# Gene to identify TF binding sites in
gene <- "NFKB2"
coordinates <- "chr10-102394110-102402524" # from UCSC site https://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=uc001kva.2&org=human

#-------------------------------------------------------------------------------
# Read in peaks
#-------------------------------------------------------------------------------

# Load database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Define list of TFs
tfs <- list.dirs("/home/ncp2306/opt/maxatac/data/models", full.names = FALSE)[2:128]

# Define function to divide res into 32bp intervals
extract_intervals <- function(row) {
  # Find number of intervals
  n <- abs(as.integer(row["end"]) - as.integer(row["start"]) - 1) / 32
  
  # Create new row for each interval
  new_row_list <- list()
  for (i in seq(n)) {
    new_row <- c(row["chromosome"],
                 as.integer(row["start"]) + 32*(i-1),
                 as.integer(row["start"]) + 32*i,
                 row["predictionScore"],
                 row["comparison"],
                 row["region_id"]
    )
    new_row <- c(new_row, paste(new_row[1], new_row[2], new_row[3], sep = "-"))
    new_row_list[[i]] <- new_row
  }
  return(do.call(rbind, new_row_list))
}

# Define frequency result outline
freqs_outline <- data.frame("Var1" = c("shared", 
                                       paste0("unique_", comparisons_oi[1]),
                                       paste0("unique_", comparisons_oi[2])))
freqs_outline$id  <- 1:nrow(freqs_outline)

# Iterate through TFs
unfinished_tfs <- c()
res_list <- list()
get_shared_unique_TFBS <- function(tf) {
  tryCatch({
    # Find files
    files <- list.files(res_dir, pattern = paste0("_", tf, "_"), full.names = TRUE)
    correct_celltype <- sapply(files, function(x) {grepl(cell_type, x, fixed = TRUE)})
    files <- files[correct_celltype]
    
    # Initialize peaks list
    peaks_ls <- list()
    
    # Read in peaks
    for (file in files) {
      # Grab comparison
      comparison <- basename(file) %>% substr(1,12)
      if (comparison %!in% comparisons_oi) {
        next
      }

      # Read in peaks
      peaks <- read.table(file, header = FALSE, sep="\t",
                          stringsAsFactors=FALSE, quote="")
      
      # Add column names
      names(peaks) <- c("chromosome", "start", "end", "predictionScore")
      
      # Isolate NFKB2
      peaks <- peaks[which(peaks$chromosome == "chr10"),]
      peaks <- peaks[which(peaks$start >= 102394110 - 3000),]
      peaks <- peaks[which(peaks$end <= 102402524 + 3000),]

      # Add comparison col
      peaks$comparison <- rep(comparison, nrow(peaks))
      
      # Append to list
      peaks_ls[[comparison]] <- peaks
    }
    
    # Check if both comparisons present
    if (length(peaks_ls) != 2) {
      unfinished_tfs <<- c(unfinished_tfs, tf)
      return(c(NA,NA,NA))
    }
    
    # Bind together
    res <- do.call(rbind, peaks_ls)
    res_list[[tf]] <<- res
    
    #-----------------------------------------------------------------------------
    # Isolate overlapping and distinct regions
    #-----------------------------------------------------------------------------

    # Extract unique chr+start+ends -> region_id
    res$region_id <- paste(res$chromosome, res$start, res$end, sep = "-")

    # Create intervals dataframe
    intervals_list <- apply(res, MARGIN = 1, FUN = extract_intervals)
    intervals <- do.call(rbind, intervals_list) %>% data.frame
    names(intervals) <- c("chromosome", "start", "end", "predictionScore", "comparison", "region_id", "interval_id")

    # Check if every interval is shared or unique
    interval_counts <- table(intervals$interval_id)
    shared_intervals <- names(interval_counts[interval_counts == 2])
    unique_intervals <- names(interval_counts[interval_counts == 1])

    # Add shared or unique info
    intervals$differential <- rep(NA, nrow(intervals))
    intervals[which(intervals$interval_id %in% shared_intervals), "differential"] <- "shared"
    intervals[which(intervals$interval_id %in% unique_intervals), "differential"] <- paste0("unique_", intervals[which(intervals$interval_id %in% unique_intervals), "comparison"])
    if (tf == "RELA") {
      rela_intervals <<- intervals
    }
    
    # Generate table
    freqs <- data.frame(table(intervals$differential))
    freqs <- merge(freqs_outline, freqs, all.x = TRUE, sort = TRUE)
    freqs <- freqs[order(freqs$id),]

    return(freqs$Freq)
  }, error = function(e) {
    message(e)
    unfinished_tfs <<- c(unfinished_tfs, tf)
    return(c(NA,NA,NA))
  })
}

# Run TFBS identification function
freqs <- lapply(tfs, get_shared_unique_TFBS)
all_regions <- do.call(rbind, res_list) %>% data.frame
saveRDS(freqs, file = paste0(output_dir, cell_type, "_ADvsHC_NFKB2_freqs.rds"))
saveRDS(all_regions, file = paste0(output_dir, cell_type, "_ADvsHC_NFKB2_regions.rds"))

#-------------------------------------------------------------------------------
# Generate visual of TF binding in NFKB2
#-------------------------------------------------------------------------------

# Merge together
freqs_df <- do.call(rbind, freqs) %>% data.frame
names(freqs_df) <- freqs_outline$Var1
freqs_df$tf <- tfs
freqs_df$total_sites <- rowSums(freqs_df[,freqs_outline$Var1], na.rm = TRUE)
freqs_df <- freqs_df[order(freqs_df$total_sites),]

# Retain only nonzero TFs 
freqs_df <- freqs_df[which(freqs_df$total_sites != 0),]

# Convert to long
freqs_long <- pivot_longer(freqs_df, cols = freqs_outline$Var1,
                           names_to = "differential", values_to = "Freq")
freqs_long$differential <- factor(freqs_long$differential, levels = freqs_outline$Var1)
freqs_long$tf <- factor(freqs_long$tf, levels = freqs_df$tf)

# Define colors
cols_diagnosis <- c("lightgray", "darkgray", "red")
cols_34vs33 <- c("lightgray", "gold", "orange")
cols_44vs33 <- c("lightgray", "gold", "red")
cols_44vs34 <- c("lightgray", "orange", "red")

# Generate plot of number share or unique
p <- ggplot(freqs_long, aes(x = Freq, y = tf, fill = differential)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.1) +
  theme_Publication_blank() +
  scale_fill_manual(values = cols_diagnosis) +
  scale_x_continuous(expand = c(0,0)) +
  # scale_x_break(c(7e5, 2e6), scales = 1, expand = FALSE) +
  labs(
    x = "# of 32bp Intervals with TFBS",
    y = "TF",
    title = paste0("Shared and unique TFBS between\n", 
                   comparisons_oi[1], " and ", 
                   comparisons_oi[2]), " in ",
                   cell_type
  ) +
  theme(axis.text.y=element_text(size=7),
        legend.position = "right")
p

# Export results
set_panel_size(p, file = paste0(output_dir, "shared_unique_TFBS_in",
                                cell_type, "_", comparisons_oi[1], "-",
                                comparisons_oi[2], "_NFKB2.pdf"),
               width = unit(4, "in"), height = unit(10, "in"))


#-------------------------------------------------------------------------------
# Grab RELA regions and plot on tracks
#-------------------------------------------------------------------------------

# Read in seurat object
s <- readRDS(noncd4_seurat_object)

# Subset for monocytes
s <- subset(s, predictedGroupNew %in% c("CD14_Mono", "CD16_Mono"))

# NFKB2 DAR
dar <- "chr10-102395377-102395877"

# define region +- 500bp around DAR
start <- 102395377 - 1000
end <- 102395877 + 1000
region <- GRanges(seqnames="chr10", ranges=paste(start, end, sep = "-"))
rela_intervals_oi <- rela_intervals[which(rela_intervals$start >= start & rela_intervals$end <= end),]

# Generate coverage plot highlighting HC binding sites
hc_rela_intervals <- rela_intervals_oi[which(rela_intervals_oi$differential == "unique_HC_Monocytes"),]
coords <- paste(hc_rela_intervals$start, hc_rela_intervals$end, sep = "-")
highlight <- GRanges(seqnames="chr10", ranges=coords)

cov_plot <- CoveragePlot(
  object = s,
  region = region,
  region.highlight = highlight,
  extend.downstream = 0,
  extend.upstream = 0,
  annotation = TRUE,
  peaks = TRUE
)

# Export results
pdf(paste0(output_dir, "RELA-binding_in_NFKB2_Monocytes_uniqueHC.pdf"))
print(cov_plot)
dev.off()

# Generate coverage plot highlighting AD binding sites
ad_rela_intervals <- rela_intervals_oi[which(rela_intervals_oi$differential == "unique_AD_Monocytes"),]
coords <- paste(ad_rela_intervals$start, ad_rela_intervals$end, sep = "-")
highlight <- GRanges(seqnames="chr10", ranges=coords)

cov_plot <- CoveragePlot(
  object = s,
  region = region,
  region.highlight = highlight,
  extend.downstream = 0,
  extend.upstream = 0,
  annotation = TRUE,
  peaks = TRUE
)

# Export results
pdf(paste0(output_dir, "RELA-binding_in_NFKB2_Monocytes_uniqueAD.pdf"))
print(cov_plot)
dev.off()

# Generate coverage plot highlighting shared binding sites
shared_rela_intervals <- rela_intervals_oi[which(rela_intervals_oi$differential == "shared"),]
shared_rela_intervals <- unique(shared_rela_intervals[,c("start", "end")])
coords <- paste(shared_rela_intervals$start, shared_rela_intervals$end, sep = "-")
highlight <- GRanges(seqnames="chr10", ranges=coords)

cov_plot <- CoveragePlot(
  object = s,
  region = region,
  region.highlight = highlight,
  extend.downstream = 0,
  extend.upstream = 0,
  annotation = TRUE,
  peaks = TRUE
)

# Export results
pdf(paste0(output_dir, "RELA-binding_in_NFKB2_Monocytes_shared.pdf"))
print(cov_plot)
dev.off()

# Generate coverage plot highlighting DAR
highlight <- GRanges(seqnames="chr10", ranges="102395377-102395877")

Idents(s) <- "Diagnosis"
cov_plot <- CoveragePlot(
  object = s,
  region = region,
  region.highlight = highlight,
  extend.downstream = 0,
  extend.upstream = 0,
  annotation = TRUE,
  peaks = TRUE
)

pdf(paste0(output_dir, "chr10-102395377-102395877-DAR_NFKB2_Monocytes.pdf"))
print(cov_plot)
dev.off()
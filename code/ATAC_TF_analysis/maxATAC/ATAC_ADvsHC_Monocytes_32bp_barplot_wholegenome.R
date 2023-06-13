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
# Summary: investigate maxatac results
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
output_dir <- "/projects/b1169/projects/AD_APOE/results_atac/maxatac/stacked_32bp_barplot_ADvsHC/out_NP_06-12-2023/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparisons of interest
comparisons_oi <- c("HC_Monocytes", "AD_Monocytes")

# Define cell type of interest
cell_type <- "Monocytes"

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
saveRDS(freqs, file = paste0(output_dir, cell_type, "_ADvsHC_freqs.rds"))
saveRDS(all_regions, file = paste0(output_dir, cell_type, "_ADvsHC_regions.rds"))

#-------------------------------------------------------------------------------
# Generate visual of everything
#-------------------------------------------------------------------------------

# Merge together
freqs <- readRDS(paste0(output_dir, cell_type, "_ADvsHC_freqs.rds"))
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
  scale_x_break(c(6e5, 1500000), scales = 0.3, expand = FALSE) +
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
pdf(paste0(output_dir, "shared_unique_TFBS_in",
           cell_type, "_", comparisons_oi[1], "-",
           comparisons_oi[2], ".pdf"),
    width = 6, height = 10)
print(p)
dev.off()


# Generate plot of number share or unique
p <- ggplot(freqs_long, aes(x = Freq, y = tf, fill = differential)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.1) +
  theme_Publication_blank() +
  scale_fill_manual(values = cols_diagnosis) +
  scale_x_continuous(expand = c(0,0)) +
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

pdf(paste0(output_dir, "shared_unique_TFBS_in",
           cell_type, "_", comparisons_oi[1], "-",
           comparisons_oi[2], "_nobreaks.pdf"),
    width = 6, height = 10)
print(p)
dev.off()
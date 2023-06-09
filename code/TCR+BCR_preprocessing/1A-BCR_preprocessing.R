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
# Written by: Abhi Ramakrishnan, modified from Natalie Piehl's scripts for TCR preprocessing
# Summary: Modify BCR barcodes and clonotype ids to match Seurat GEX formatting
# Create .csv containing filtered contig annotations for all samples in AD-APOE project
# Create .csv containing paired (both alpha and beta sequences) clonotype annotations for all samples in AD-APOE project
# Modified from original script by Emma Tapp from Stanford University
# ------------------------------------------------------------------------------

## Load libraries and define variables
# Load packages
library(tidyverse)
library(purrr)
library(data.table)
library("Seurat")
library("plyr")
library(tidyr)
# Organize directories
bcr_dir <- "/path/to/cellranger/bcr/outputs/"
output_dir <- "/path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
#################################################################################################################
# Merge and modify contig annotations
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}
# Make a list of cellranger bcr output dirs for each sample
out_dirs <- unlist(lapply(list.dirs.depth.n(bcr_dir, 2), function(x) grep("*/outs", x, value = TRUE)))

# Create list of paths to filtered_contig_annotations.csv files
filtered_contig_paths <- list.files(path = out_dirs, pattern = "filtered_contig_annotations.csv", full.names = TRUE, recursive = TRUE)

# Initialize list to contain contig dataframes
filtered_contig_dfs <- list()

# Generate list of dataframes of filtered_contig_annotations.csv files
for (path in filtered_contig_paths) {
  # Read in csv
  tmp <- read.csv(path)

  # Isolate sample ID from path
  id <- unlist(strsplit(path, "/")) %>%
    tail(3) %>%
    pluck(1)

  # Create ID  column
  tmp$id <- rep(id, dim(tmp)[1])

  # Add dataframe to list
  filtered_contig_dfs[[id]] <- tmp
}

# Merge all dataframes into one
filtered_contig_merged <- rbindlist(filtered_contig_dfs)
head(filtered_contig_merged)

# Add sample ID to barcode and raw_clonotype_id
filtered_contig_merged$barcode <- paste0(filtered_contig_merged$id,"_",filtered_contig_merged$barcode)
filtered_contig_merged$raw_clonotype_id <- paste0(filtered_contig_merged$id,"_",filtered_contig_merged$raw_clonotype_id)

# Remove cells with no consensus clonotype
filtered_contig_merged <- filtered_contig_merged[which(filtered_contig_merged$raw_consensus_id != 'None'),]

# Write out merged contigs
write.csv(filtered_contig_merged, file = paste0(output_dir, "/contigs_merged.csv"), row.names = FALSE)

#################################################################################################################
## Merge and modify clonotype annotations

# Create list of paths to Clonotype<ID>.csv files
clonotype_paths <- list.files(path = out_dirs, pattern = "clonotypes.csv", full.names = TRUE, recursive = TRUE)

# Initialize list to cointain contig dataframes
clonotype_dfs <- list()

# Generate list of dataframes of filtered_contig_annotations.csv files
for (path in clonotype_paths) {
  # Read in csv
  tmp <- read.csv(path)

  # Isolate sample name from path
  id <- unlist(strsplit(path, "/")) %>%
    tail(3) %>%
    pluck(1)

  # Create id column
  tmp$id <- id

  # Add dataframe to list
  clonotype_dfs[[id]] <- tmp
}

# Merge all dataframes into one
clonotype_merged <- rbindlist(clonotype_dfs)
head(clonotype_merged)

# Break apart CDR3 column into separate consensus sequences
clonotype_merged <- separate(data = clonotype_merged, col = cdr3s_aa, into = c("seq1", "seq2", "seq3", "seq4"), sep = ";")
# Remove rows that had more than 4 sequences
clonotype_merged <- clonotype_merged[-c(2627,2628,2649,2660),]
# Create empty columns for IGH and IGK/L sequences
clonotype_merged$igh_cdr3s <- rep(NA, dim(clonotype_merged)[1])
clonotype_merged$igk_cdr3s <- rep(NA, dim(clonotype_merged)[1])
clonotype_merged$igl_cdr3s <- rep(NA, dim(clonotype_merged)[1])

# Remove unneccesary columns
clonotype_merged$cdr3s_nt <- NULL

for (row_num in seq(nrow(clonotype_merged))) {
  # Create IGH column containing <igh_consensus_1>;<igh_consensus_2>
  row <- clonotype_merged[row_num,]
  igh_seqs <- grep("^IGH:", row, value = TRUE) %>%
    str_remove_all("IGH:") %>%
    paste(collapse = ";")

  clonotype_merged$igh_cdr3s[ row_num ] <- igh_seqs

  # Create IGK column containing <igk_consensus_1>;<igk_consensus_2>
  igk_seqs <- grep("^IGK:", row, value = TRUE) %>%
    str_remove_all("IGK:") %>%
    paste(collapse = ";")

  clonotype_merged$igk_cdr3s[ row_num ] <- igk_seqs

  # Create IGL column containing <igl_consensus_1>;<igl_consensus_2>
  igl_seqs <- grep("^IGL:", row, value = TRUE) %>%
    str_remove_all("IGL:") %>%
    paste(collapse = ";")

  clonotype_merged$igl_cdr3s[ row_num ] <- igl_seqs
}

clonotype_merged <- subset(clonotype_merged, select=-c(seq1, seq2, seq3, seq4))
head(clonotype_merged)

# Isolate paired clonotypes
paired_clonotypes_merged <- clonotype_merged[which(clonotype_merged$igh_cdr3s != "" & (clonotype_merged$igk_cdr3s != "" | clonotype_merged$igl_cdr3s != "" )),]

# Add sample ID to clonotype_id
paired_clonotypes_merged$clonotype_id <- paste0(paired_clonotypes_merged$id,"_",paired_clonotypes_merged$clonotype_id)

# Write out merged paired clonotypes
write.csv(paired_clonotypes_merged, paste0(output_dir, "/paired_clonotypes_merged.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Add bcr data (lines 108-145 of merging with Seurat object script)
#-------------------------------------------------------------------------------
# Load in latest seurat object (which contains tcr metadata)
scRNA <- "path/to/object/"
scRNA <- load(scRNA)

# Load in bcr information
contigs_merged <- read.csv(paste0(output_dir, "/contigs_merged.csv"))
bcrs_paired <- read.csv(paste0(output_dir, "/paired_clonotypes_merged.csv"))

## Now add in bcr information to meta data
# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
bcr <- contigs_merged[!duplicated(contigs_merged$barcode), ]
# Add letter G to front of each bcr barcode to match seurat object barcodes
bcr$barcode <- paste0("G", bcr$barcode)
# Only keep the barcode and clonotype columns
# We'll get additional clonotype info from the clonotype table.
bcr <- bcr[,c("barcode", "raw_clonotype_id")]
names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"

# Slap the AA sequences onto our original table by clonotype_id.
bcr <- merge(bcr, bcrs_paired[, c("clonotype_id", "igh_cdr3s", "igk_cdr3s","igl_cdr3s", "frequency")])
head(bcr)
# Reorder so barcodes are first column and set them as rownames.
bcr <- bcr[, c(2,1,3,4,5,6)]
rownames(bcr) <- bcr[,1]
bcr[,1] <- NULL
colnames(bcr) <- c("B_clonotype_id", "igh_cdr3s", "igk_cdr3s", "igl_cdr3s", "B_frequency")
print("Head of bcr")
head(bcr)
# Add to object
head(s@meta.data)
s <- AddMetaData(object=s,
                 metadata=bcr)
print("Head of metadata")
head(s@meta.data)
table(s@meta.data$B_frequency)

#------------------------------------------------------------------------------
# Filter bcrs (modified from lines 39-73 of Natalie's tcr_cleanup-main script)
#------------------------------------------------------------------------------
# Remove bcr annotations from Non B cells
table(s$predicted.celltype.l2)
b_cell_clusters <- c(
  "B naive",
  "B intermediate",
  "B memory",
  "Plasmablast"
)
s@meta.data$B_clonotype_id[
  which(s@meta.data$predicted.celltype.l2 %!in% b_cell_clusters)] <- NA
s@meta.data$igh_cdr3s[
  which(s@meta.data$predicted.celltype.l2 %!in% b_cell_clusters)] <- NA
s@meta.data$igk_cdr3s[
  which(s@meta.data$predicted.celltype.l2 %!in% b_cell_clusters)] <- NA
s@meta.data$igl_cdr3s[
  which(s@meta.data$predicted.celltype.l2 %!in% b_cell_clusters)] <- NA
s@meta.data$B_frequency[
  which(s@meta.data$predicted.celltype.l2 %!in% b_cell_clusters)] <- NA

print(table(s[[c("predicted.celltype.l2", "B_frequency")]]))

# Calculate filtered frequency
filtered_freq <- data.frame(table(s[["B_clonotype_id"]]))
names(filtered_freq) <- c("B_clonotype_id", "frequency_filtered")
filtered_freq_merged <- plyr::join(s[["B_clonotype_id"]], filtered_freq)

# If in same order, merge
if (all.equal(filtered_freq_merged$B_clonotype_id, s[["B_clonotype_id"]]$B_clonotype_id)) {
  s@meta.data$B_frequency_filtered <- filtered_freq_merged$frequency_filtered
} else {
  stop("Unable to merge.")
}
print(table(s[[c("predicted.celltype.l2", "B_frequency_filtered")]]))
# ------------------------------------------------------------------------------
# Add clonal info to the Seurat object's metadata. (Modified lines 167-205 of Natalie's preprocessing-Seurat script)
# ------------------------------------------------------------------------------
# Add clonal info to metadata
s@meta.data$clonal_B <- "NA"
s@meta.data$clonal_B[
  s@meta.data$B_frequency_filtered > 1] <- "C"
s@meta.data$clonal_B[
  s@meta.data$B_frequency_filtered == 1] <- "NC"
table(s@meta.data$clonal_B)
print(table(s[[c("predicted.celltype.l2", "clonal_B")]]))

# Save modified seurat object
save(s, file = paste0(output_dir, "s_VDJclean"))
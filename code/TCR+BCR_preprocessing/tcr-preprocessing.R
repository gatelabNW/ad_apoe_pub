# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 08-30-2022
# Written by: Natalie Piehl, modified for AD_APOE project by Abhi Ramakrishnan
# Summary: Modify barcodes and clonotype ids to match Seurat GEX formatting
# Create .csv containing filtered contig annotations for all samples in CSF aging project
# Create .csv containing paired (both alpha and beta sequences) clonotype annotations for all samples in CSF aging project
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
tcr_dir <- "/path/to/cellranger/tcr/outputs/"
output_dir <- "/path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
#################################################################################################################
## Merge and modify contig annotations
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}
# Make a list of cellranger tcr output dirs for each sample
out_dirs <- unlist(lapply(list.dirs.depth.n(tcr_dir, 2), function(x) grep("*/outs", x, value = TRUE)))

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
# Remove a row that had more than 4 sequences
clonotype_merged <- clonotype_merged[-12679,]
# Create empty columns for TRA and TRB sequences
clonotype_merged$tra_cdr3s <- rep(NA, dim(clonotype_merged)[1])
clonotype_merged$trb_cdr3s <- rep(NA, dim(clonotype_merged)[1])

# Remove unneccesary columns
clonotype_merged$cdr3s_nt <- NULL
clonotype_merged$inkt_evidence <- NULL
clonotype_merged$mait_evidence <- NULL

for (row_num in seq(nrow(clonotype_merged))) {
  # Create TRA column containing <tra_consensus_1>;<tra_consensus_2>
  row <- clonotype_merged[row_num,]
  tra_seqs <- grep("^TRA:", row, value = TRUE) %>%
    str_remove_all("TRA:") %>%
    paste(collapse = ";")

  # clonotype_merged[ row_num, "tra_cdr3s"] <- tra_seqs
  clonotype_merged$tra_cdr3s[ row_num ] <- tra_seqs

  # Create TRB column containing <trb_consensus_1>;<trb_consensus_2>
  trb_seqs <- grep("^TRB:", row, value = TRUE) %>%
    str_remove_all("TRB:") %>%
    paste(collapse = ";")

  clonotype_merged$trb_cdr3s[ row_num ] <- trb_seqs
}

clonotype_merged <- subset(clonotype_merged, select=-c(seq1, seq2, seq3, seq4))
head(clonotype_merged)

# Isolate paired clonotypes
paired_clonotypes_merged <- clonotype_merged[which(clonotype_merged$trb_cdr3s != "" & clonotype_merged$tra_cdr3s != "" & clonotype_merged$tra_cdr3s != 0 ),]

# Add sample ID to clonotype_id
paired_clonotypes_merged$clonotype_id <- paste0(paired_clonotypes_merged$id,"_",paired_clonotypes_merged$clonotype_id)

# Write out merged paired clonotypes
write.csv(paired_clonotypes_merged, paste0(output_dir, "/paired_clonotypes_merged.csv"), row.names = FALSE)

#-------------------------------------------------------------------------------
# Add TCR data (lines 108-145 of merging with Seurat object script)
#-------------------------------------------------------------------------------
# Load seurat object with RNA data
scRNA <- "path/to/object/"
scRNA <- load(scRNA)
# Add cell ID to metadata
ID <- sub("(.*?_.*?)_.*", "\\1", rownames(s@meta.data))
s <- AddMetaData(object=s,
                 metadata=data.frame(ID=ID,
                                     row.names=rownames(s@meta.data)))
# Load in TCR information
contigs_merged <- read.csv(paste0(output_dir, "/contigs_merged.csv"))
TCRs_paired <- read.csv(paste0(output_dir, "/paired_clonotypes_merged.csv"))

## Now Add in TCR information to meta data
# Subsets so only the first line of each barcode is kept,
# as each entry for given barcode will have same clonotype.
tcr <- contigs_merged[!duplicated(contigs_merged$barcode), ]
# Add letter G to front of each tcr barcode to match seurat object barcodes
tcr$barcode <- paste0("G", tcr$barcode)
# Only keep the barcode and clonotype columns
# We'll get additional clonotype info from the clonotype table.
tcr <- tcr[,c("barcode", "raw_clonotype_id")]
names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

# Slap the AA sequences onto our original table by clonotype_id.
tcr <- merge(tcr, TCRs_paired[, c("clonotype_id", "trb_cdr3s", "tra_cdr3s", "frequency")])

# Reorder so barcodes are first column and set them as rownames.
tcr <- tcr[, c(2,1,3,4, 5)]
rownames(tcr) <- tcr[,1]
tcr[,1] <- NULL
colnames(tcr) <- c("T_clonotype_id", "trb_cdr3s", "tra_cdr3s", "T_frequency")
print("Head of TCR")
head(tcr)

s <- AddMetaData(object=s,
                 metadata=tcr)
print("Head of metadata")
head(s@meta.data)
table(s@meta.data$T_frequency)

#------------------------------------------------------------------------------
# Filter TCRs (modified from lines 39-73 of Natalie's tcr_cleanup-main script)
#------------------------------------------------------------------------------
## Remove TCR annotations from Non T cells
# Define T cell types
t_cell_clusters <-c("CD4 CTL",
                    "CD4 Naive",
                    "CD4 Proliferating",
                    "CD4 TCM",
                    "CD4 TEM",
                    "Treg",
                    "CD8 Naive",
                    "CD8 Proliferating",
                    "CD8 TCM",
                    "CD8 TEM",
                    "MAIT",
                    "dnT",
                    "gdT")

s@meta.data$T_clonotype_id[
  which(s@meta.data$predicted.celltype.l2 %!in% t_cell_clusters)] <- NA
s@meta.data$trb_cdr3s[
  which(s@meta.data$predicted.celltype.l2 %!in% t_cell_clusters)] <- NA
s@meta.data$tra_cdr3s[
  which(s@meta.data$predicted.celltype.l2 %!in% t_cell_clusters)] <- NA
s@meta.data$T_frequency[
  which(s@meta.data$predicted.celltype.l2 %!in% t_cell_clusters)] <- NA

print(table(s[[c("predicted.celltype.l2", "T_frequency")]]))

# Calculate filtered frequency
filtered_freq <- data.frame(table(s[["T_clonotype_id"]]))
names(filtered_freq) <- c("T_clonotype_id", "frequency_filtered")
filtered_freq_merged <- plyr::join(s[["T_clonotype_id"]], filtered_freq)

# If in same order, merge
if (all.equal(filtered_freq_merged$T_clonotype_id, s[["T_clonotype_id"]]$T_clonotype_id)) {
  s@meta.data$T_frequency_filtered <- filtered_freq_merged$frequency_filtered
} else {
  stop("Unable to merge.")
}
print(table(s[[c("predicted.celltype.l2", "T_frequency_filtered")]]))
# ------------------------------------------------------------------------------
# Add clonal info to the Seurat object's metadata. (Modified lines 167-205 of Natalie's preprocessing-Seurat script)
# ------------------------------------------------------------------------------
# Add clonal info to metadata
s@meta.data$clonal_T <- "NA"
s@meta.data$clonal_T[
  s@meta.data$T_frequency_filtered > 1] <- "C"
s@meta.data$clonal_T[
  s@meta.data$T_frequency_filtered == 1] <- "NC"
table(s@meta.data$clonal_T)
print(table(s[[c("predicted.celltype.l2", "clonal_T")]]))

## Add clonality with diagnosis
# Initialize diagnosis.clonal_T column in metadata
s@meta.data$diagnosis.clonal_T <- "NA"

# AD
tmp <- s@meta.data[which(s@meta.data$Diagnosis=="Alzheimers Disease"),]
AD.C <- tmp[which(tmp$clonal == "C"),]
AD.NC <- tmp[which(tmp$clonal == "NC"),]
s@meta.data$diagnosis.clonal_T[
  which(rownames(s@meta.data) %in% rownames(AD.C))] <- "AD.C"
s@meta.data$diagnosis.clonal_T[
  which(rownames(s@meta.data) %in% rownames(AD.NC))] <- "AD.NC"

# HC
tmp <- s@meta.data[which(s@meta.data$Diagnosis=="Healthy Control"),]
HC.C <- tmp[which(tmp$clonal == "C"),]
HC.NC <- tmp[which(tmp$clonal == "NC"),]
s@meta.data$diagnosis.clonal_T[
  which(rownames(s@meta.data) %in% rownames(HC.C))] <- "HC.C"
s@meta.data$diagnosis.clonal_T[
  which(rownames(s@meta.data) %in% rownames(HC.NC))] <- "HC.NC"

table(s@meta.data$diagnosis.clonal_T)

#------------------------------------------------------------------------------
# Generate normalized frequency
#------------------------------------------------------------------------------
# Calculate TCR counts per sample
tcr_count <- s[[c("predicted.celltype.l2", "T_clonotype_id", "Age", "T_frequency_filtered")]] %>%
  dplyr::filter(!is.na(T_frequency_filtered)) %>%
  dplyr::distinct(T_clonotype_id, .keep_all = TRUE) %>%
  dplyr::group_by(predicted.celltype.l2) %>%
  dplyr::summarize(tcr_count = sum(T_frequency_filtered), across(Age)) %>%
  dplyr::distinct()
head(tcr_count)

# Calculate median TCR counts
median_tcr_count <- median(tcr_count$tcr_count)
median_tcr_count

# Initialize column
s@meta.data$normalized_frequency_T <- "NA"

# Normalize frequency for each sample
for (sample in tcr_count$predicted.celltype.l2) {
  sample_tcr_count <- tcr_count %>% dplyr::filter(predicted.celltype.l2 == sample) %>% dplyr::pull(tcr_count)
  s@meta.data$normalized_frequency_T[
    which(s@meta.data$predicted.celltype.l2 == sample) ] <- s@meta.data$T_frequency_filtered[
      which(s@meta.data$predicted.celltype.l2 == sample) ] * (median_tcr_count / sample_tcr_count)
}

# Check output
head(s[[c("T_frequency", "T_frequency_filtered", "normalized_frequency_T")]])

# Save modified seurat object
save(s, file = paste0(output_dir, "s_tcrclean"))

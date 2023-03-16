# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-15-2022
# Written by: Natalie Piehl
# Summary: Generate coverage plot for region of interest
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("GenomicRanges")
})

# read in both CD4 and non-CD4 objects
cd4 <- readRDS("/path/to/object/")
non_cd4 <- readRDS("/path/to/object/")

# List groups of RNA fine celltypes
#-------------------------------------------------------------------------------
rna_B <- c(
  "B_naive",
  "B_intermediate",
  "B_memory",
  "Plasmablast"
)
rna_CD4_T <- c(
  "CD4_CTL",
  "CD4_Proliferating",
  "CD4_Naive",
  "CD4_TCM",
  "CD4_TEM",
  "Treg"
)
rna_CD8_T <- c(
  "CD8_Naive",
  "CD8_Proliferating",
  "CD8_TCM",
  "CD8_TEM",
  "MAIT"
)
rna_NK <- c(
  "NK",
  "NK_CD56bright",
  "NK_Proliferating"
)
rna_Mono <- c(
  "CD14_Mono",
  "CD16_Mono"
)
rna_dc <- c(
  "cDC2",
  "pDC",
  "ASDC"
)

rna_other_T <- c(
  "ILC",
  "gdT",
  "dnT"
)
rna_other <- c(
  "Platelet",
  "Eryth",
  "HSPC"
)
#-------------------------------------------------------------------------------

# set output directory
output_dir <- "/path/to/output/directory/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Subset for celltype of interest
s <- subset(non_cd4, predictedGroupNew == rna_Mono)

# Make APOE + Diagnosis column
s@meta.data$Diagnosis_APOE <- paste(s@meta.data$Diagnosis, s@meta.data$APOE_genotype)

# Define levels
s@meta.data$Diagnosis_APOE <- factor(s@meta.data$Diagnosis_APOE,
                                     levels = c("Healthy Control E3/E3",
                                                "Alzheimers Disease E3/E3",
                                                "Healthy Control E3/E4",
                                                "Alzheimers Disease E3/E4",
                                                "Healthy Control E4/E4",
                                                "Alzheimers Disease E4/E4"))
s@meta.data$Diagnosis <- factor(s@meta.data$Diagnosis,
                                levels = c("Healthy Control",
                                           "Alzheimers Disease"))
s@meta.data$APOE_genotype <- factor(s@meta.data$APOE_genotype,
                                    levels = c("E3/E3",
                                               "E3/E4",
                                               "E4/E4"))

# Downsample so each group has same number of cells
barcodes <- s[["Diagnosis_APOE"]]
barcodes$barcode <- rownames(barcodes)
freq <- table(barcodes$Diagnosis_APOE) %>% as.data.frame()
sample_num <- min(freq$Freq)
barcodes_sampled <- barcodes %>%
  group_by(Diagnosis_APOE) %>%
  slice_sample(n = sample_num, replace = FALSE)
s <- subset(s, cells = barcodes_sampled$barcode)

#-------------------------------------------------------------------------------
# Specify gene of interest
gene <- "IL1B"
#-------------------------------------------------------------------------------

# 1.) Plot IL1B_1 distal region
# Degine range of interest
range_oi <- GRanges(seqnames="chr2", ranges= "112811894-112812394")

# Set identity
Idents(s) <- "Diagnosis"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = GRanges(seqnames = "chr2", ranges="112811894-112812394"),
  region.highlight = range_oi,
  extend.downstream = 500,
  extend.upstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  # ymax = 100
)

pdf(paste0(output_dir, "Monocyte_", gene,"_1_distal_region", "_diagnosis_covplot.pdf"))
print(cov_plot)
dev.off()

# 2.) Plot IL1B_2 promoter region
# Define range of interest
range_oi <- GRanges(seqnames="chr2", ranges= "112836690-112837190")

# Set identity
Idents(s) <- "Diagnosis"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = gene,
  region.highlight = range_oi,
  extend.downstream = 500,
  extend.upstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  # ymax = 100
)

pdf(paste0(output_dir, "Monocyte_", gene,"_2_promoter_region", "_diagnosis_covplot.pdf"))
print(cov_plot)
dev.off()

# 3.) Plot IL1B_3 distal region
# Define range of interest
range_oi <- GRanges(seqnames="chr2", ranges= "112839288-112839788")

# Set identity
Idents(s) <- "Diagnosis"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = range_oi,
  region.highlight = range_oi,
  extend.downstream = 500,
  extend.upstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  # ymax = 100
)

pdf(paste0(output_dir, "Monocyte_", gene,"_3_distal_region", "_diagnosis_covplot.pdf"))
print(cov_plot)
dev.off()

#-------------------------------------------------------------------------------
# Specify gene of interest
gene <- "CCL4L2"
#-------------------------------------------------------------------------------

# 4.) Plot CCL4L2_1 promoter AND CCL4L2_2 exonic region
# Define range of interest
range_oi <- GRanges(seqnames="chr17", ranges= c("36210694-36211194", "36211850-36212350" ))

# Set identity
Idents(s) <- "Diagnosis"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = gene,
  region.highlight = range_oi,
  extend.downstream = 500,
  extend.upstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  # ymax = 100
)

pdf(paste0(output_dir, "Monocyte_", gene,"_1_promoter+", gene, "_2_exonic_regions", "_diagnosis_covplot.pdf"))
print(cov_plot)
dev.off()

#-------------------------------------------------------------------------------
# Specify gene of interest
gene <- "CCL3"
#-------------------------------------------------------------------------------

# 5.) Plot CCL3_3 promoter region
# Define range of interest
range_oi <- GRanges(seqnames="chr17", ranges= "36089944-36090444")

# Set identity
Idents(s) <- "Diagnosis"

# Generate coverage plot
cov_plot <- CoveragePlot(
  object = s,
  region = gene,
  region.highlight = range_oi,
  extend.downstream = 500,
  extend.upstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  # ymax = 100
)

pdf(paste0(output_dir, "Monocyte_", gene,"_3_promoter_region", "_diagnosis_covplot.pdf"))
print(cov_plot)
dev.off()


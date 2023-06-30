# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 02-16-2023
# Written by: Natalie Piehl
# Summary: Run Cicero
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("monocle3")
  library("cicero")
  library("Seurat")
  library("Signac")
  library("SeuratWrappers")
})

# Organize inputs
ranges_path <- "/path/to/ranges/"
input_base_dir <- "/path/to/cicero/results"
output_dir <- "/path/to/output_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Set up
#-------------------------------------------------------------------------------

# Read in ranges
ranges <- readRDS(ranges_path)

# Load in STAR gene data
ref_path <- "/path/to/ensembl/Homo_sapiens.GRCh38.109.chr.gtf.gz"
gene_anno <- rtracklayer::readGFF(ref_path)

# rename some columns to match plotting requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#-------------------------------------------------------------------------------
# Generate connection plots
#-------------------------------------------------------------------------------

# Define cell type
cell_type <- "CD8+_T_Cells"

# Load in data
diagnosis <- "Alzheimers_Disease"
input_dir <- paste0(input_base_dir, diagnosis, "/")
conns_ad <- readRDS(paste0(input_dir, cell_type, "_conns.rds"))
ccans_ad <- readRDS(paste0(input_dir, cell_type, "_ccans.rds"))

diagnosis <- "Healthy_Control"
input_dir <- paste0(input_base_dir, diagnosis, "/")
conns_hc <- readRDS(paste0(input_dir, cell_type, "_conns.rds"))
ccans_hc <- readRDS(paste0(input_dir, cell_type, "_ccans.rds"))

# CXCR3
pdf(file = paste0(output_dir, "CXCR3_chrX-71623759-71624259_", cell_type, "_ADvsHC_connections.pdf"))
plot_connections(conns_hc,
                 comparison_track = conns_ad,
                 chr = "chrX",
                 viewpoint = "chrX_71623759_71624259",
                 connection_ymax = 0.3,
                 comparison_ymax = 0.3,
                 collapseTranscripts = "longest",
                 comparison_connection_color = "red",
                 connection_color = "gray",
                 alpha_by_coaccess = TRUE,
                 minbp = 71600000,
                 maxbp = 71630000,
                 gene_model = gene_anno)
dev.off()

# CXCR3
pdf(file = paste0(output_dir, "CXCR3_chrX-71619128-71619628_", cell_type, "_ADvsHC_connections.pdf"))
plot_connections(conns_hc,
                 comparison_track = conns_ad,
                 chr = "chrX",
                 viewpoint = "chrX_71619128_71619628",
                 connection_ymax = 0.3,
                 comparison_ymax = 0.3,
                 collapseTranscripts = "longest",
                 comparison_connection_color = "red",
                 connection_color = "gray",
                 alpha_by_coaccess = TRUE,
                 minbp = 71610000,
                 maxbp = 71630000,
                 gene_model = gene_anno)
dev.off()

# Define cell type
cell_type <- "Monocytes"

# Load in data
diagnosis <- "Alzheimers_Disease"
input_dir <- paste0(input_base_dir, diagnosis, "/")
conns_ad <- readRDS(paste0(input_dir, cell_type, "_conns.rds"))
ccans_ad <- readRDS(paste0(input_dir, cell_type, "_ccans.rds"))

diagnosis <- "Healthy_Control"
input_dir <- paste0(input_base_dir, diagnosis, "/")
conns_hc <- readRDS(paste0(input_dir, cell_type, "_conns.rds"))
ccans_hc <- readRDS(paste0(input_dir, cell_type, "_ccans.rds"))


# ABCA1 in Monocytes
pdf(file = paste0(output_dir, "ABCA1_promoter_chr9-104927943-104928443_", cell_type, "_ADvsHC_connections.pdf"))
plot_connections(conns_hc,
                 comparison_track = conns_ad,
                 chr = "chr9",
                 viewpoint = "chr9_104927943_104928443",
                 connection_ymax = 0.3,
                 comparison_ymax = 0.3,
                 collapseTranscripts = "longest",
                 comparison_connection_color = "red",
                 connection_color = "gray",
                 alpha_by_coaccess = TRUE,
                 minbp = 104900000,
                 maxbp = 104935000,
                 gene_model = gene_anno)
dev.off()

# ABCA1 in Monocytes
pdf(file = paste0(output_dir, "ABCA1_intron_chr9-104927943-104928443_", cell_type, "_ADvsHC_connections.pdf"))
plot_connections(conns_hc,
                 comparison_track = conns_ad,
                 chr = "chr9",
                 viewpoint = "chr9_104909664_104910164",
                 connection_ymax = 0.3,
                 comparison_ymax = 0.3,
                 collapseTranscripts = "longest",
                 comparison_connection_color = "red",
                 connection_color = "gray",
                 alpha_by_coaccess = TRUE,
                 minbp = 104900000,
                 maxbp = 104935000,
                 gene_model = gene_anno)
dev.off()

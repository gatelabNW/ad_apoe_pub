# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-23-2023
# Written by: Natalie Piehl
# Summary: Run Footprinting
#
#-------------------------------------------------------------------------------
# Install packages

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("optparse")
})

# Generate list of expect arguments 
option_list = list(
  make_option(c("-c", "--celltype"), action="store", type="character", default=NULL, 
              help="cell type to analyze"),
    make_option(c("-m", "--motif"), action="store", type="character", default=NULL, 
              help="motif to analyze")
); 
 
# Parse these arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Handle null arguments
if (is.null(opt$celltype)){
  print_help(opt_parser)
  stop("A cell type must be supplied", call.=FALSE)
} else if (is.null(opt$motif)){
  print_help(opt_parser)
  stop("A motif must be supplied", call.=FALSE)
}

# Save arguments to variable
cell_type <- opt$celltype
print(cell_type)
motif <- opt$motif
print(motif)

# Organize inputs
output_dir <- paste0("/projects/b1169/projects/AD_APOE/results_atac/footprinting/batch/out_NP_02-09-2023/", cell_type, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Run footprinting

# Load in seurat object
if (cell_type == "CD4+_T_Cells") {
  s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds"))
} else {
  s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"))
}
s

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

# Subset for celltype of interest
s <- subset(s, broad_celltype == cell_type)
s

# Add APOE_Diagnosis column
s@meta.data$APOE_Diagnosis <- paste(s@meta.data$APOE_genotype, s@meta.data$Diagnosis, sep = "_")

# Set identify to APOE x Diagnosis
ident <- "Diagnosis"
s <- SetIdent(s, value = ident)

# Add footprinting
s <- Footprint(
  object = s,
  motif.name = c(motif),
  in.peaks = TRUE,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Generate plot
p <- PlotFootprint(s, features = motif)
pdf(paste0(output_dir, cell_type, "_", motif, "_", ident, "_footprint.pdf"))
print(p)
dev.off()

# Save footprint output
saveRDS(s@assays$peaks@positionEnrichment, 
        file = paste0(output_dir, cell_type, "_", motif, "_positionEnrichment.rds"))
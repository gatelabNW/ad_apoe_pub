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
  library("optparse")
})

# Generate list of expect arguments 
option_list = list(
  make_option(c("-d", "--diagnosis"), action="store", type="character", default=NULL, 
              help="diagnosis to analyze"),
  make_option(c("-g", "--geno"), action="store", type="character", default=NULL, 
              help="genotype to analyze"),
  make_option(c("-t", "--celltype"), action="store", type="character", default=NULL, 
              help="cell type to analyze")
); 

# Parse these arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Handle null arguments
if (is.null(opt$diagnosis)){
  print_help(opt_parser)
  stop("A diagnosis must be supplied", call.=FALSE)
} else if (is.null(opt$celltype)){
  print_help(opt_parser)
  stop("A cell type must be supplied", call.=FALSE)
} else if (is.null(opt$geno)){
  print_help(opt_parser)
  stop("A genotype must be supplied", call.=FALSE)
}

# Save arguments to variable
diagnosis <- opt$diagnosis
geno <- opt$geno
cell_type <- opt$celltype
print(diagnosis)
print(geno)
print(cell_type)

# Organize inputs
output_base_dir <- "/projects/b1169/projects/AD_APOE/results_atac/cicero/batch/out_NP_02-20-2023/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Set up
#-------------------------------------------------------------------------------

# Generate new output dir
if (geno == "all") {
  output_dir <- paste0(output_base_dir, 
                       gsub(" ", "_", diagnosis), "/")
} else {
  output_dir <- paste0(output_base_dir, 
                       gsub(" ", "_", diagnosis), "_", 
                       gsub("/", "", geno), "/")
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load in seurat object
print("Loading and formatting seurat object...")
if (cell_type == "CD4+_T_Cells") {
  s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds"))
  
  # Read in UMAP and add
  # s <- readRDS("/projects/b1169/projects/AD_APOE/results/archr/manual_outs/to_signac/2022_11_29_NP/cd4_s_withUMAP.rds")
  # umap <- s@reductions$umap
  # saveRDS(umap, file = paste0(output_base_dir, "cd4_umap.rds"))
  umap <- readRDS(paste0(output_base_dir, "cd4_umap.rds"))
  s@reductions$umap <- umap
} else {
  s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"))
  
  # Read in UMAP and add
  # s <- readRDS("/projects/b1169/projects/AD_APOE/results/archr/manual_outs/to_signac/2022_11_29_NP/noncd4_s_withUMAP.rds")
  # umap <- s@reductions$umap
  # saveRDS(umap, file = paste0(output_base_dir, "noncd4_umap.rds"))
  umap <- readRDS(paste0(output_base_dir, "noncd4_umap.rds"))
  s@reductions$umap <- umap
}

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

# Subset for celltype and diagnosis
s <- subset(s, broad_celltype == cell_type)
if (diagnosis != "all") {
  s <- subset(s, Diagnosis == diagnosis)
}
if (geno != "all") {
  s <- subset(s, APOE_genotype == geno)
}

# Check for any peaks with no accessibility in celltype/diagnosis 
count_matrix <- GetAssayData(s, slot='counts', assay='peaks')
count_matrix <- count_matrix[rowSums(count_matrix) > 0,]
s <- subset(s, features = rownames(count_matrix))
s

# # Make CDS with SeuratWrappers
# # 6.8 GB
# # contains raw and norm counts
cds <- SeuratWrappers::as.cell_data_set(x = s)

# Make CDS with Monocle3
# 3.5 GB
# Contains only raw counts
# print("Generating cds...")
# genes <- data.frame(as.character(rownames(count_matrix)))
# rownames(genes) <- rownames(count_matrix)
# genes <- as.data.frame(cbind(genes,genes))
# colnames(genes) <- c("GeneSymbol", "gene_short_name")
# cds <- new_cell_data_set(
#   count_matrix,
#   cell_metadata=s@meta.data,
#   gene_metadata=genes
# )

# Make Cicero CDS
cicero <- make_cicero_cds(cds,
                          reduced_coordinates = reducedDims(cds)$UMAP,
                          k = 50)

# Export Cicero CDS object
saveRDS(cicero, file = paste0(output_dir, cell_type, "_cicero-cds.rds"))

# Get the chromosome sizes from the Seurat object
genome <- seqlengths(s)

# Convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

#-------------------------------------------------------------------------------
# Run cicero
#-------------------------------------------------------------------------------

# Run cicero
print("Calculating connections...")
conns <- run_cicero(cicero,
                    genomic_coords = genome.df,
                    sample_num = 100)

print("Calculating CCANs...")
# Generate CCANs
ccans <- generate_ccans(conns)

# Export results
saveRDS(conns, file = paste0(output_dir, cell_type, "_conns.rds"))
saveRDS(ccans, file = paste0(output_dir, cell_type, "_ccans.rds"))
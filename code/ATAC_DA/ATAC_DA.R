# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 01-30-2023
# Written by: Natalie Piehl
# Summary: Run DA
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
  library("optparse")
})

# Generate list of expect arguments 
option_list = list(
  make_option(c("-c", "--comparison"), action="store", type="character", default=NULL, 
              help="comparison to make"),
  make_option(c("-t", "--celltype"), action="store", type="character", default=NULL, 
              help="cell type to analyze")
); 

# Parse these arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Handle null arguments
if (is.null(opt$comparison)){
  print_help(opt_parser)
  stop("A comparison must be supplied", call.=FALSE)
} else if (is.null(opt$celltype)){
  print_help(opt_parser)
  stop("A cell type must be supplied", call.=FALSE)
}

# Save arguments to variable
comparison <- opt$comparison
cell_type <- opt$celltype
print(comparison)
print(cell_type)

# Organize inputs
ranges_path <- "/path/to/granges(seurat_object)"
output_base_dir <- "/path/to/output/folder"
output_dir <- paste0(output_base_dir, comparison, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

#-------------------------------------------------------------------------------
# Format data
#-------------------------------------------------------------------------------

# Load in seurat object
if (cell_type == "CD4+_T_Cells") {
  s <- readRDS("/path/to/cd4_seurat_object")
} else {
  s <- readRDS("/path/to/noncd4_seurat_object")
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

# Read in ranges
ranges <- readRDS(ranges_path)

# Subset for cell type
if (cell_type != "all") {
  print(paste0("Subsetting for ", cell_type))
  s <- subset(s, broad_celltype == cell_type)
}

# Set ident and subset based on comparison
if ((comparison != "ADvsHC") & grepl("ADvsHC", comparison)) {
  # Obtain genotype
  genotype <- gsub("ADvsHC_", "", comparison)
  genotype_split <- strsplit(genotype, "") %>% unlist
  genotype_formatted <- paste0("E", genotype_split[1],
                               "/", "E", genotype_split[2])
  
  # Subset for genotype
  print(paste0("Subsetting for ", genotype_formatted))
  s <- subset(s, APOE_genotype == genotype_formatted)
  
  # Set Idents
  print(paste0("Setting ident to Diagnosis: AD over HC"))
  s <- SetIdent(s, value = "Diagnosis")
  ident_1 <- "Alzheimers Disease"
  ident_2 <- "Healthy Control"
} else if (comparison == "ADvsHC") {
  # Set Idents
  print(paste0("Setting ident to Diagnosis: AD over HC"))
  s <- SetIdent(s, value = "Diagnosis")
  ident_1 <- "Alzheimers Disease"
  ident_2 <- "Healthy Control"
} else {
  # Obtain genotypes to compare
  genotype_1 <- str_sub(comparison, start = 1, end = 2)
  genotype_2 <- str_sub(comparison, start = 5, end = 6)
  genotype_1_split <- strsplit(genotype_1, "") %>% unlist
  genotype_2_split <- strsplit(genotype_2, "") %>% unlist
  ident_1 <- paste0("E", genotype_1_split[1],
                    "/", "E", genotype_1_split[2])
  ident_2 <- paste0("E", genotype_2_split[1],
                    "/", "E", genotype_2_split[2])
  
  # Obtain diagnosis
  diagnosis <- str_sub(comparison, start = -2)
  if (diagnosis == "AD") {
    diagnosis_formatted <- "Alzheimers Disease"
  } else {
    diagnosis_formatted <- "Healthy Control"
  }
  
  # Subset for diagnosis
  print(paste0("Subsetting for ", diagnosis_formatted))
  s <- subset(s, Diagnosis == diagnosis_formatted)
  
  # Set Idents
  print(paste0("Setting ident to APOE Genotype: ", ident_1, " over ", ident_2))
  s <- SetIdent(s, value = "APOE_genotype")
}

# Verify correct samples retained
print("Retained cell type(s), diagnosis(es), and APOE genotype(s)...")
print(unique(s[["predictedGroupNew"]])[,1])
print(unique(s[["Diagnosis"]])[,1])
print(unique(s[["APOE_genotype"]])[,1])

#-------------------------------------------------------------------------------
# Run DA
#-------------------------------------------------------------------------------

# Define output csv path
output_csv_path <- paste0(output_dir, cell_type, "_", comparison, "_dars.csv")

# Check if analysis has been run...
if (!file.exists(output_csv_path)) {
  tryCatch({
    # Print status update
    print("Running DA now...")
    
    # Specify covariates
    if (comparison == "ADvsHC") {
      covars <- c('nCount_peaks', 'Sex', 'APOE_genotype')
    } else {
      covars <- c('nCount_peaks', 'Sex')
    }
    print(paste0("Using the following covariates: ", paste(covars, collapse = ", ")))
    
    # Run DA
    da_peaks <- FindMarkers(
      object = s,
      ident.1 = ident_1,
      ident.2 = ident_2,
      only.pos = FALSE,
      logfc.threshold = -Inf,
      test.use = 'LR',
      min.pct = 0.05,
      latent.vars = covars
    )

    # Run BH procedure
    da_peaks$BH <- p.adjust(da_peaks$p_val, method = "BH")

    # Annotate with genes
    da_peaks$sitename <- rownames(da_peaks)
    da_peaks <- merge(da_peaks, ranges, all.x = TRUE, all.y = FALSE)

    # Annotate with numbered gene region column
    da_peaks$gene_region_index <- rep(NA, nrow(da_peaks))
    gene_list <- list()
    for (i in seq(nrow(da_peaks))) {
      if (da_peaks[i, "nearestGene"] %!in% names(gene_list)) {
        gene_list[[da_peaks[i, "nearestGene"]]] <- 1
      } else {
        gene_list[[da_peaks[i, "nearestGene"]]] <- gene_list[[da_peaks[i, "nearestGene"]]] + 1
      }
      da_peaks[i, "gene_region_index"] <- paste(da_peaks[i, "nearestGene"],
                                                gene_list[[da_peaks[i, "nearestGene"]]],
                                                sep = "_")
    }
    
    # Export DARs
    write.csv(da_peaks, file = output_csv_path)

    # Change format to match volcano script
    da_peaks <- da_peaks[which(da_peaks$gene_region_index != "NA_"),]
    rownames(da_peaks) <- da_peaks$gene_region_index

    # Generate volcano plot
    volcano_plot(da_peaks, file = paste0(output_dir, cell_type, "_", comparison, "_volcano.pdf"),
                 lfc.thresh = lfc.thresh, padj.thresh = padj.thresh,
                 title = paste0(comparison, " in ", cell_type))

  }, error = function(e) {
    message(e)
    NULL
  })
} else {
  print(paste0(comparison, " DA has already been run on ", cell_type))
}
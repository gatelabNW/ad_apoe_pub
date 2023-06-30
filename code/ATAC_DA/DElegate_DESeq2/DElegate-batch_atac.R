# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-05-2023
# Written by: Natalie Piehl
# Summary: Run DA with DElegate
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Increase max size
options(future.globals.maxSize=1048576000000)

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Signac")
  library("DElegate")
  library("optparse")
})

# Generate list of expect arguments 
option_list = list(
  make_option(c("-c", "--comparison"), action="store", type="character", default=NULL,
              help="comparison to make")
); 

# Parse these arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Handle null arguments
if (is.null(opt$comparison)){
  print_help(opt_parser)
  stop("A comparison must be supplied", call.=FALSE)
}

# Save arguments to variable
comparison <- opt$comparison
print(comparison)

# Organize inputs
celltype_colors_path <- "/path/to/celltype_color_map.csv"
output_dir <- paste0("/path/to/output/dir/", comparison, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

#-------------------------------------------------------------------------------
# Set up (nonCD4 first)
#-------------------------------------------------------------------------------

# Load in seurat object
s_atac <- readRDS(paste0("/path/to/noncd4/seurat/object"))

# Create new objecta
meta <- s_atac[[]]
meta$RNA <- NULL
s <- CreateSeuratObject(assay = "RNA", 
                        counts = s_atac@assays$peaks@counts,
                        meta.data = meta)
rm(s_atac)

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

# Make E4 carrier identity class
s@meta.data$E4_carrier <- mapvalues(s@meta.data$APOE_genotype,
                                    from = c("E3/E3", "E3/E4", "E4/E4"),
                                    to = c("noncarrier", "carrier", "carrier"))

# Set ident and subset based on comparison
if (grepl("ADvsHC", comparison)) {
  # Obtain genotype
  genotype <- gsub("ADvsHC_", "", comparison)
  genotype_split <- strsplit(genotype, "") %>% unlist
  genotype_formatted <- paste0("E", genotype_split[1],
                               "/", "E", genotype_split[2])
  
  if (comparison == "ADvsHC_4carrier") {
    # Subset for genotype
    print(paste0("Subsetting for E4 carriers"))
    s <- subset(s, E4_carrier == "carrier") 
  } else if (comparison != "ADvsHC") {
    # Subset for genotype
    print(paste0("Subsetting for ", genotype_formatted))
    s <- subset(s, APOE_genotype == genotype_formatted) 
  }
  
  # Define idents
  ident_1 <- "Alzheimers Disease"
  ident_2 <- "Healthy Control"
  group <- "Diagnosis"
} else if (grepl("carrier", comparison)) {
  if (comparison != "4carriervs33") {
    # Obtain diagnosis
    diagnosis <- str_sub(comparison, start = -2)
    if (diagnosis == "AD") {
      diagnosis_formatted <- "Alzheimers Disease"
    } else {
      diagnosis_formatted <- "Healthy Control"
    }
    
    # Subset for diagnosis and APOE
    print(paste0("Subsetting for ", diagnosis_formatted))
    s <- subset(s, Diagnosis == diagnosis_formatted)
  }
  
  # Define idents
  ident_1 <- "carrier"
  ident_2 <- "noncarrier"
  group <- "E4_carrier"
} else {
  # Obtain genotypes to compare
  genotype_1 <- str_sub(comparison, start = 1, end = 2)
  genotype_2 <- str_sub(comparison, start = 5, end = 6)
  genotype_1_split <- strsplit(genotype_1, "") %>% unlist
  genotype_2_split <- strsplit(genotype_2, "") %>% unlist
  ident_1 <- paste0("E", genotype_1_split[1],
                    "/E", genotype_1_split[2])
  ident_2 <- paste0("E", genotype_2_split[1],
                    "/E", genotype_2_split[2])
  group <- "APOE_genotype"
  
  if (grepl("AD", comparison) | grepl("HC", comparison)) {
    # Obtain diagnosis
    diagnosis <- str_sub(comparison, start = -2)
    if (diagnosis == "AD") {
      diagnosis_formatted <- "Alzheimers Disease"
    } else {
      diagnosis_formatted <- "Healthy Control"
    }
    
    # Subset for diagnosis and APOE
    print(paste0("Subsetting for ", diagnosis_formatted, " in ", ident_1, " and ", ident_2))
    s <- subset(s, Diagnosis == diagnosis_formatted & APOE_genotype %in% c(ident_1, ident_2))
  } else {
    # Subset for APOE
    print(paste0("Subsetting for ", ident_1, " and ", ident_2))
    s <- subset(s, APOE_genotype %in% c(ident_1, ident_2))
  }
}

# Verify correct samples retained
print("Retained diagnosis(es), and APOE genotype(s)...")
print(unique(s[["Diagnosis"]])[,1])
print(unique(s[["APOE_genotype"]])[,1])
print(group)
print(ident_1)
print(ident_2)

#-------------------------------------------------------------------------------
# Run DE
#-------------------------------------------------------------------------------

run_de <- function(cell_type) {
  tryCatch({
    print(cell_type)
    cell_type_label <- gsub("/", "", cell_type)
    cell_type_label <- gsub(" ", "_", cell_type_label)
    
    # Extract data of interest
    s_celltype <- subset(s, broad_celltype == cell_type)
    
    # Isolate genes expressed in >5% of cells
    data <-  GetAssayData(object = s_celltype, assay = "RNA", slot = "data")
    data <- data[ which((rowSums(data > 0) / ncol(data)) > .05), ]
    s_celltype <- subset(s_celltype, features = rownames(data))
    
    # Run DE
    de_results <- findDE(object = s_celltype,
                         replicate_column = "Sample",
                         group_column = group,
                         method = "deseq",
                         compare = c(ident_1, ident_2))
    
    # Write out results
    write.csv(de_results, 
              file = paste0(output_dir, cell_type_label, "_", comparison, "_deseq_dars.csv"))
  }, error = function(e) {
    message(e)
  })
}

# Run DE on all celltypes
cell_types <- unique(s[["broad_celltype"]])[,1]
lapply(cell_types, run_de)

#-------------------------------------------------------------------------------
# Set up (CD4 next)
#-------------------------------------------------------------------------------

# Load in seurat object
s_atac <- readRDS(paste0("/path/to/cd4/seurat/object"))

# Create new objecta
meta <- s_atac[[]]
meta$RNA <- NULL
s <- CreateSeuratObject(assay = "RNA", 
                        counts = s_atac@assays$peaks@counts,
                        meta.data = meta)
rm(s_atac)

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

# Make E4 carrier identity class
s@meta.data$E4_carrier <- mapvalues(s@meta.data$APOE_genotype,
                                    from = c("E3/E3", "E3/E4", "E4/E4"),
                                    to = c("noncarrier", "carrier", "carrier"))

# Set ident and subset based on comparison
if (grepl("ADvsHC", comparison)) {
  # Obtain genotype
  genotype <- gsub("ADvsHC_", "", comparison)
  genotype_split <- strsplit(genotype, "") %>% unlist
  genotype_formatted <- paste0("E", genotype_split[1],
                               "/", "E", genotype_split[2])
  
  if (comparison == "ADvsHC_4carrier") {
    # Subset for genotype
    print(paste0("Subsetting for E4 carriers"))
    s <- subset(s, E4_carrier == "carrier") 
  } else if (comparison != "ADvsHC") {
    # Subset for genotype
    print(paste0("Subsetting for ", genotype_formatted))
    s <- subset(s, APOE_genotype == genotype_formatted) 
  }
  
  # Define idents
  ident_1 <- "Alzheimers Disease"
  ident_2 <- "Healthy Control"
  group <- "Diagnosis"
} else if (grepl("carrier", comparison)) {
  if (comparison != "4carriervs33") {
    # Obtain diagnosis
    diagnosis <- str_sub(comparison, start = -2)
    if (diagnosis == "AD") {
      diagnosis_formatted <- "Alzheimers Disease"
    } else {
      diagnosis_formatted <- "Healthy Control"
    }
    
    # Subset for diagnosis and APOE
    print(paste0("Subsetting for ", diagnosis_formatted))
    s <- subset(s, Diagnosis == diagnosis_formatted)
  }
  
  # Define idents
  ident_1 <- "carrier"
  ident_2 <- "noncarrier"
  group <- "E4_carrier"
} else {
  # Obtain genotypes to compare
  genotype_1 <- str_sub(comparison, start = 1, end = 2)
  genotype_2 <- str_sub(comparison, start = 5, end = 6)
  genotype_1_split <- strsplit(genotype_1, "") %>% unlist
  genotype_2_split <- strsplit(genotype_2, "") %>% unlist
  ident_1 <- paste0("E", genotype_1_split[1],
                    "/E", genotype_1_split[2])
  ident_2 <- paste0("E", genotype_2_split[1],
                    "/E", genotype_2_split[2])
  group <- "APOE_genotype"
  
  if (grepl("AD", comparison) | grepl("HC", comparison)) {
    # Obtain diagnosis
    diagnosis <- str_sub(comparison, start = -2)
    if (diagnosis == "AD") {
      diagnosis_formatted <- "Alzheimers Disease"
    } else {
      diagnosis_formatted <- "Healthy Control"
    }
    
    # Subset for diagnosis and APOE
    print(paste0("Subsetting for ", diagnosis_formatted, " in ", ident_1, " and ", ident_2))
    s <- subset(s, Diagnosis == diagnosis_formatted & APOE_genotype %in% c(ident_1, ident_2))
  } else {
    # Subset for APOE
    print(paste0("Subsetting for ", ident_1, " and ", ident_2))
    s <- subset(s, APOE_genotype %in% c(ident_1, ident_2))
  }
}

# Verify correct samples retained
print("Retained diagnosis(es), and APOE genotype(s)...")
print(unique(s[["Diagnosis"]])[,1])
print(unique(s[["APOE_genotype"]])[,1])
print(group)
print(ident_1)
print(ident_2)

# Run DE on all celltypes
cell_types <- unique(s[["broad_celltype"]])[,1]
lapply(cell_types, run_de)

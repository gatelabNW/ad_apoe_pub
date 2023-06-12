# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-13-2022
# Written by: Morabito, Natalie Piehl
# Summary: Follow Morabito code to correlate RNA expression and ATAC
# accessibility
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
  library("pbapply")
})

# Generate list of expect arguments 
option_list = list(
  make_option(c("-a", "--atac_celltype"), action="store", type="character", default=NULL, 
              help="ATAC cell type to analyze"),
  make_option(c("-r", "--rna_celltype"), action="store", type="character", default=NULL, 
              help="RNA cell type to analyze"),
  make_option(c("-c", "--coacc_cutoff"), action="store", type="double", default=NULL, 
              help="coacc cutoff to use"),
  make_option(c("-d", "--diagnosis"), action="store", type="character", default=NULL, 
              help="diagnosis to use"),
  make_option(c("-p", "--apoe"), action="store", type="character", default=NULL, 
              help="apoe to use")
); 

# Parse these arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Handle null arguments
if (is.null(opt$atac_celltype)){
  print_help(opt_parser)
  stop("An ATAC cell type must be supplied", call.=FALSE)
} else if (is.null(opt$rna_celltype)){
  print_help(opt_parser)
  stop("An RNA cell type must be supplied", call.=FALSE)
} else if (is.null(opt$coacc_cutoff)){
  print_help(opt_parser)
  stop("A coaccessibility cutoff must be supplied", call.=FALSE)
}

# Save arguments to variable
atac_cell_type <- opt$atac_celltype
rna_cell_type <- opt$rna_celltype
coacc_cutoff <- opt$coacc_cutoff
diagnosis <- opt$diagnosis
apoe <- opt$apoe
message(paste0("Comparing acc of ", atac_cell_type, " with exp of ", rna_cell_type,
             " using coacc cutoff of ", coacc_cutoff, " in ", diagnosis, " ", apoe, " samples"))

# Organize inputs
rna_seurat_object <-"/projects/b1169/projects/AD_APOE/results/seurat/supervised_clustering/out_NP_09-07-2022/s_sup_clustering"
cicero_dir <- "/projects/b1169/projects/AD_APOE/results_atac/cicero/batch/out_NP_02-20-2023/"
output_dir <- paste0("/projects/b1169/projects/AD_APOE/results_atac/expression_accessibility_correlation/batch/out_NP_04-03-2023/", atac_cell_type, "-", rna_cell_type, "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Prep
#-------------------------------------------------------------------------------

# Load RNA seurat object and ATAC proj + seurat object
message("---------- Loading seurat objects... ----------")
load(rna_seurat_object)
if (atac_cell_type == "CD4+_T_Cells") {
  atac_s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/cd4_s_TFIDF.rds"))
} else {
  atac_s <- readRDS(paste0("/projects/b1169/projects/AD_APOE/results_atac/conversion/TFIDF_normalization/out_NP_02-06-2023/noncd4_s_TFIDF.rds"))
}

# Add broad cell types
atac_s@meta.data$broad_celltype <- mapvalues(atac_s@meta.data$predictedGroupNew,
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

# Subset RNA and ATAC seurat objects by celltype and reserve only shared samples
message("---------- Subsetting for shared samples... ----------")
shared_samples <- intersect(unique(s[["orig.ident"]])[,1],
                            sapply(unique(atac_s[["Sample"]])[,1],
                                   function(x) {gsub("A", "G", x)}))
s <- subset(s, orig.ident %in% shared_samples)
atac_s <- subset(atac_s, Sample %in% sapply(shared_samples,
                                            function(x) {gsub("G", "A", x)}))

# Set idents to sampleIDs:
Idents(s) <- "orig.ident"
Idents(atac_s) <- "Sample"
DefaultAssay(s) <- "RNA"
DefaultAssay(atac_s) <- "peaks"

#-------------------------------------------------------------------------------
# Define functions
#-------------------------------------------------------------------------------

# Format seurat objects
format_seurat <- function(cell_type, diagnosis, apoe, type) {
  # Subset seurat objects for cell type
  if (type == "RNA") {
    subset <- subset(s, predicted.celltype.l2 == gsub("_", " ", cell_type))
  } else if (type == "ATAC") {
    subset <- subset(atac_s, broad_celltype == cell_type)
  } else {
    message("Invalid type given. Please provide either RNA or ATAC")
  }
  
  # Subset for diagnosis and apoe if needed
  if (!is.null(diagnosis)) {
    subset <- subset(subset, Diagnosis == gsub("_", " ", diagnosis))
  } 
  if (!is.null(apoe)) {
    subset <- subset(subset, APOE_genotype == apoe)
  }
  return(subset)
}

# Format cicero objects
format_conns <- function(cell_type, diagnosis, apoe,
                         rna_subset, atac_subset) {
  # Load conns for conditions
  if (is.null(diagnosis)) {
    if (is.null(apoe)) {
      conns_path <- paste0(cicero_dir, "all/", cell_type, "_conns.rds")
    } else {
      conns_path <- paste0(cicero_dir, "all_", gsub("/", "", apoe), "/", cell_type, "_conns.rds")
    }
  } else {
    if (is.null(apoe)) {
      conns_path <- paste0(cicero_dir, diagnosis, "/", cell_type, "_conns.rds")
    } else {
      conns_path <- paste0(cicero_dir, diagnosis, "_", gsub("/", "", apoe), "/", cell_type, "_conns.rds")
    }
  }
  conns <- readRDS(conns_path)
  
  # Isolate granges object data for all linked peaks
  ranges <- granges(atac_subset) %>% as.data.frame
  ranges$site_name <- paste(ranges$seqnames, ranges$start, ranges$end, sep = "-")
  tmp1 <- ranges[na.omit(match(conns$Peak1, ranges$site_name)),
                         c('peakType', 'nearestGene', 'site_name')]
  tmp2 <- ranges[na.omit(match(conns$Peak2, ranges$site_name)),
                         c('peakType', 'nearestGene', 'site_name')]
  
  # add columns to conns
  conns$Peak1_type <- tmp1$peakType
  conns$Peak1_nearestGene <- tmp1$nearestGene
  conns$Peak2_type <- tmp2$peakType
  conns$Peak2_nearestGene <- tmp2$nearestGene
  
  # Remove conns where nearest gene is NA
  conns <- conns[!is.na(conns$Peak1_nearestGene),]

  # get only links that are connected to a gene promoter region
  conns <- subset(conns, Peak1_type == 'Promoter')

  # split conns by target gene:
  conns_list <- group_split(conns, Peak1_nearestGene)
  names(conns_list) <- sort(unique(conns$Peak1_nearestGene))
  
  return(conns_list)
}

run_corr <- function(gene, conns_list, avg_exp, avg_acc) {
  # Subset for gene of interest and co-accessibility
  conns <- conns_list[[gene]]
  
  # Remove any conns between promoters of the same gene
  conns <- conns[!(conns$Peak2 %in% conns$Peak1),]
  
  # Subset by coaccessibility score
  conns <- subset(conns, coaccess >= coacc_cutoff)

  # Skip this gene if there are no co-accessible connections
  if(nrow(conns) == 0){return(data.frame())}
  
  # Retain only connections where accessibility of second peak is nonzero
  # in at least 5 samples
  conns <- conns[which(conns$Peak2 %in% rownames(avg_acc)),]

  # get average exp and acc for this gene and peaks that are co-accessible
  avg_acc_gene <- avg_acc[as.character(conns$Peak2),]
  avg_exp_gene <- avg_exp[gene,] %>% as.vector
  
  # Skip this gene if no peak present
  if(is.null(dim(avg_acc_gene))){
    print(paste0("No peak present for ", gene))
    return(data.frame())}

  # correlation between expression and accessibility:
  cor_mat <- apply(avg_acc_gene, 1, function(x){
    correlation <- cor.test(as.numeric(avg_exp_gene),
                            as.numeric(x),
                            method='pearson')
    data.frame("pval"=correlation$p.value,
               "pcc"=correlation$estimate,
               "avg_exp"=paste(sapply(avg_exp_gene, function(x) {round(x, 3)}), collapse = ", "),
               "avg_acc"=paste(sapply(x, function(y) {round(y, 3)}), collapse = ", "))
  })

  # collapse individual correlation dfs
  cor_df <- Reduce(rbind, cor_mat)

  # add correlation stats to conns
  conns$pcc <- cor_df$pcc
  conns$pval <- cor_df$pval
  conns$avg_exp <- cor_df$avg_exp
  conns$avg_acc <- cor_df$avg_acc
  
  return(conns)
}

#-------------------------------------------------------------------------------
# Execute
#-------------------------------------------------------------------------------

# Define main function
main <- function(atac_cell_type, rna_cell_type, diagnosis, apoe) {
  # Print status
  print(paste(atac_cell_type, rna_cell_type, diagnosis))
  
  # Format seurat objects
  message("---------- Formatting seurat objects... ----------")
  rna_subset <- format_seurat(cell_type = rna_cell_type,
                              diagnosis = diagnosis,
                              apoe = apoe,
                              type = "RNA")
  atac_subset <- format_seurat(cell_type = atac_cell_type,
                               diagnosis = diagnosis,
                               apoe = apoe,
                               type = "ATAC")
  
  # Subset for shared samples within subtype
  message("---------- Subsetting for shared samples within celltype... ----------")
  shared_samples <- intersect(unique(rna_subset[["orig.ident"]])[,1],
                              sapply(unique(atac_subset[["Sample"]])[,1],
                                     function(x) {gsub("A", "G", x)}))
  rna_subset <- subset(rna_subset, orig.ident %in% shared_samples)
  atac_subset <- subset(atac_subset, Sample %in% sapply(shared_samples,
                                              function(x) {gsub("G", "A", x)}))
  
  # Compute average exp/acc
  avg_acc <- AverageExpression(atac_subset)$peaks
  avg_exp <- AverageExpression(rna_subset)$RNA
  
  # Remove genes without expression in at least 5 samples
  avg_exp <- avg_exp[which(rowSums(avg_exp != 0) >= 5),]
  
  # Remove regions without accessibility in at least 5 samples
  avg_acc <- avg_acc[which(rowSums(avg_acc != 0) >= 5),]

  # Ensure samples are in same order
  atac_sample_order <- sapply(colnames(avg_acc), function(x) {gsub("A", "G", x)})
  avg_exp <- avg_exp[,atac_sample_order]
  
  # Format cicero data
  message("---------- Formatting cicero data... ----------")
  conns_list <- format_conns(cell_type = atac_cell_type,
                             diagnosis = diagnosis,
                             apoe = apoe,
                             rna_subset = rna_subset,
                             atac_subset = atac_subset)
  
  # Extract genes to analyze
  rna_genes <- sapply(rownames(avg_exp), function(x) {gsub("\\.","-",x)})
  genes <- names(conns_list)[names(conns_list) %in% rna_genes]
  
  # Run correlation
  message(paste0("---------- Running correlation on ", length(genes), " genes... ----------"))
  corr_list <- pblapply(genes, run_corr,
                        conns_list = conns_list,
                        avg_acc = avg_acc,
                        avg_exp = avg_exp)

  # combine into one df and remove incomplete entries:
  message("---------- Merging results... ----------")
  df <- do.call(rbind, corr_list)
  df <- df[!is.na(df$pcc),]
  
  # compute corrected p vals
  message("---------- Calculating adjusted p values... ----------")
  df$BH <- p.adjust(df$pval, method='BH')

  # write df to file:
  saveRDS(df, file=paste0(output_dir, diagnosis, "_", gsub("/", "", apoe), "_", 
                          atac_cell_type, "-", rna_cell_type,
                          '_peak_gene_correlation.rds'))
  message("---------- Finished! ----------")
}

# Run all
main(atac_cell_type = atac_cell_type, rna_cell_type = rna_cell_type,
     diagnosis = diagnosis, apoe = apoe)

#-------------------------------------------------------------------------------
# Analyze
#-------------------------------------------------------------------------------

# Read in results
df <- readRDS(paste0(output_dir, diagnosis, "_", gsub("/", "", apoe), "_", 
                     atac_cell_type, "-", rna_cell_type,
                     '_peak_gene_correlation.rds'))

# Subset for 95th percentile of PCC and adjusted p val < .01
pcc_thresh <- quantile(df$pcc, 0.95)
padj_thresh <- .01
df <- df[which(df$pcc >= pcc_thresh & df$BH <= padj_thresh),]

# Change commas to spaces
df$avg_exp <- sapply(df$avg_exp, function(x) {gsub(",", " ", x)})
df$avg_acc <- sapply(df$avg_acc, function(x) {gsub(",", " ", x)})

# Export sig correlations
write.table(df,
            file=paste0(output_dir, diagnosis, "_", gsub("/", "", apoe), "_", 
                        atac_cell_type, "-", rna_cell_type,
                        '_peak_gene_sig_correlation.csv'),
            sep=',', quote=FALSE, row.names=FALSE)

# Get number of gl-cCREs and cCRE linked genes
num_cres <- length(unique(df$Peak2))
num_genes <- length(unique(df$Peak1_nearestGene))
print(paste0("gl-cCREs: ", num_cres))
print(paste0("cCRE linked genes: ", num_genes))

# Find numbers of gl-cCREs per cCRE linked gene
freq <- table(df$Peak1_nearestGene) %>% as.data.frame

# Make barplot of distribution
p <- ggplot(freq, aes(x = Freq)) +
  geom_histogram(binwidth=1, colour="black", fill = "darkolivegreen1") +
  scale_y_continuous(expand = c(0.03,0.03)) +
  scale_x_continuous(expand = c(0.03,0.03), breaks = seq(max(freq$Freq))) +
  theme_Publication_blank() +
  labs(title = paste0(diagnosis, "\n", atac_cell_type, "-", rna_cell_type),
       x = "Number of gl-cCREs",
       y = "Number of cCRE linked genes") +
  theme(text = element_text(size=20))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, gsub(" ", "_", diagnosis), "_", gsub("/", "", apoe), "_", 
                                atac_cell_type, "-", rna_cell_type,
                                '_num_sig_glCRE_distribution.pdf'))

# Find numbers of peak types
peak_types <- table(df$Peak2_type) %>% as.data.frame
print(peak_types)
# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 09-19-2022
# Written by: Natalie Piehl
# Summary: Run DE on Diagnosis
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("doMC")
  library("UpSetR")
})

# Organize inputs
celltype_colors_path <- "/projects/b1169/nat/als/resources/metadata/cluster-metadata.csv"
seurat_object <- "/projects/b1169/projects/AD_APOE/results/seurat/supervised_clustering/out_NP_09-07-2022/s_sup_clustering"
output_dir <- "/projects/b1169/projects/AD_APOE/results/de/diagnosis/out_NP_09-28-2022_covarSex+APOE/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# Set core number for parallel model fitting
registerDoMC(cores = 2)

# Load Seurat object
load(seurat_object)

#-------------------------------------------------------------------------------
# Run DE on AD vs HC

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set Ident to Diagnosis
s <- SetIdent(s, value = "Diagnosis")
print(table(s[["Diagnosis"]]))

run_de <- function(cell_type) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)

  # Find DEGs b/w MCI/AD and HC
  degs <-FindMarkers(object = subset(s, predicted.celltype.l2 == cell_type),
                     ident.1 = "Alzheimers Disease",
                     ident.2 = "Healthy Control",
                     latent.vars = c("Sex", "APOE_genotype"),
                     test.use = "MAST",
                     logfc.threshold = -Inf,
                     min.pct = 0.1,
                     assay = "RNA"
  )

  # Remove ribosomal, mitochondrial, and HLA genes
  if (length(grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs))) != 0) {
    degs <- degs[-grep(pattern = "^RPS|^RPL|^MT-|^HLA-", x = rownames(degs)),] 
  }

  # Run Benjamini-Hochberg adjustment
  degs$BH <- p.adjust(degs$p_val, method = "BH")
  print(head(degs))

  # Write out results
  write.csv(degs, paste0(output_dir, cell_type_label, "_AD_vs_HC_degs.csv"))

  # Create volcano plot
  volcano_plot(degs, title = paste0("AD vs HC in ", cell_type),
               file = paste0(output_dir, cell_type_label, "_AD_vs_HC_volcano.pdf"),
               padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
  print(paste0(cell_type, " is done!"))
}

# # Run DE on all celltypes
cell_types <- unique(s[["predicted.celltype.l2"]])[,1]
mclapply(cell_types, run_de, mc.cores = 2)
# lapply(cell_types, run_de)

#------------------------------------------------------------------------------
# Create Upset plot
sig_genes_ls <- list()

# Create list with sig genes for each cell type
for (cell_type in cell_types) {
  print(cell_type)
  cell_type_label <- gsub("/", "", cell_type)
  cell_type_label <- gsub(" ", "_", cell_type_label)
  
  # Load in degs
  tryCatch({
    degs <- read.csv(paste0(output_dir, cell_type_label, "_AD_vs_HC_degs.csv"))
    
    # Identify sig genes
    sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
    print(head(sig_genes))
    
    # Add sig genes to list
    sig_genes_ls[[cell_type]] <- sig_genes$X
  }, error = function(e) {
    NULL
  })
  
}

sig_genes_ls
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get colors
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"
    
# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]
celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "/upset_celltype.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()
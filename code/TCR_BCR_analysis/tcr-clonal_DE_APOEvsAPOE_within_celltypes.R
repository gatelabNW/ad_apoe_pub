# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 12-14-2022
# Written by: Abhi Ramakrishnan
# Summary: Perform differential expression analysis of clonal T-cells by 
# diagnosis and genotype within celltypes for AD-APOE project
# Sex covariate
# ------------------------------------------------------------------------------
# Create output directory
output_dir <- "path/to/output/directory/"
ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# Source Natalie's helper functions
source("path/to/helper_functions.R")

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("UpSetR")
  library(ggplot2)
})
# Load seurat object
s <- "path/to/object/"
load(s)

# Run standard normalization
DefaultAssay(object = s) <- "RNA"
s <- NormalizeData(s, verbose = FALSE)

# Set identity used for DE
s <- SetIdent(s, value = "APOE_genotype")

# Subset Seurat object by clonality
object_clonal <- subset(s, clonal_T== "C")

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

# List diagnoses
diag_list <- c("Healthy Control", "Alzheimers Disease")

# List celltypes
# CD4 proliferating, CD8 Naive, CD8 proliferating, and dnT no longer a part of this list due to low cell number
# excluding all same celltypes for each comparison for consistency
t_cell_clusters <-c("CD4 CTL",
                    "CD4 Naive",
                    "CD4 TCM",
                    "CD4 TEM",
                    "Treg",
                    # "CD8 Naive",
                    # "CD8 Proliferating",
                    "CD8 TCM",
                    "CD8 TEM",
                    "MAIT",
                    "gdT")
# ------------------------------------------------------------------------------
# 1. Run DE for each genotype comparison within each diagnosis subset
# ------------------------------------------------------------------------------
# Run DE for APOE 3/4 vs 3/3
for (diag in diag_list) {
  print(diag)
  # Subset Seurat object by diagnosis
  object_diag <- subset(object_clonal, Diagnosis== diag)
  for (cluster in t_cell_clusters) {
    x <- FindMarkers(object = subset(object_diag, predicted.celltype.l2==cluster),
                     ident.1 = "E3/E4",
                     ident.2= "E3/E3",
                     logfc.threshold = -Inf,
                     test.use = "MAST",
                     min.pct = 0.1,
                     assay= 'RNA',
                     latent.vars = 'Sex'
    )
    # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
    if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
      x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
    }
    # Run Benjamini-Hochberg adjustment
    x$BH <- p.adjust(x$p_val, method = "BH")
    print(head(x))

    # Write out results
    diag <- gsub(" ", "_", diag)
    cluster_name <- gsub(" ", "_", cluster)
    write.csv(x, paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe34_vs_33_degs.csv"))

    # Create volcano plot
    volcano_plot(x, title = paste0("APOE E3/E4 vs E3/E3 in ", diag, " ", cluster_name),
                 file = paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe34_vs_33_volcano.pdf"),
                 padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
    print(paste0(diag,"_", cluster_name, " apoe 3/4 vs 3/3 is done!"))
  }

}

# Run DE for APOE 4/4 vs 3/4

for (diag in diag_list) {
  print(diag)
  # Subset Seurat object by diagnosis
  object_diag <- subset(object_clonal, Diagnosis== diag)
  for (cluster in t_cell_clusters) {
    x <- FindMarkers(object = subset(object_diag, predicted.celltype.l2==cluster),
                     ident.1 = "E4/E4",
                     ident.2= "E3/E4",
                     logfc.threshold = -Inf,
                     test.use = "MAST",
                     min.pct = 0.1,
                     assay= 'RNA',
                     latent.vars = 'Sex'
    )
    # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
    if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
      x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
    }
    # Run Benjamini-Hochberg adjustment
    x$BH <- p.adjust(x$p_val, method = "BH")
    print(head(x))

    # Write out results
    diag <- gsub(" ", "_", diag)
    cluster_name <- gsub(" ", "_", cluster)
    write.csv(x, paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_34_degs.csv"))

    # Create volcano plot
    volcano_plot(x, title = paste0("APOE E4/E4 vs E3/E4 in ", diag, " ", cluster_name),
                 file = paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_34_volcano.pdf"),
                 padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
    print(paste0(diag,"_", cluster_name, " apoe 4/4 vs 3/4 is done!"))
  }

}

# Run DE for APOE 4/4 vs 3/3

for (diag in diag_list) {
  print(diag)
  # Subset Seurat object by diagnosis
  object_diag <- subset(object_clonal, Diagnosis== diag)
  for (cluster in t_cell_clusters) {
    x <- FindMarkers(object = subset(object_diag, predicted.celltype.l2==cluster),
                     ident.1 = "E4/E4",
                     ident.2= "E3/E3",
                     logfc.threshold = -Inf,
                     test.use = "MAST",
                     min.pct = 0.1,
                     assay= 'RNA',
                     latent.vars = 'Sex'
    )
    # Remove ribosomal, mitochondrial (leaving HLA genes present for now)
    if (length(grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x))) != 0) {
      x <- x[-grep(pattern = "^RPS|^RPL|^MT-", x = rownames(x)),]
    }
    # Run Benjamini-Hochberg adjustment
    x$BH <- p.adjust(x$p_val, method = "BH")
    print(head(x))

    # Write out results
    diag <- gsub(" ", "_", diag)
    cluster_name <- gsub(" ", "_", cluster)
    write.csv(x, paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_33_degs.csv"))

    # Create volcano plot
    volcano_plot(x, title = paste0("APOE E4/E4 vs E3/E3 in ", diag, " ", cluster_name),
                 file = paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_33_volcano.pdf"),
                 padj.thresh = padj.thresh, lfc.thresh = lfc.thresh)
    print(paste0(diag,"_", cluster_name, " apoe 4/4 vs 3/3 is done!"))
  }

}

# # ------------------------------------------------------------------------------
# # 2. Make upset plots for each diagnosis group to quantify num DEGs per genotype comparison
# # ------------------------------------------------------------------------------

# A.) UPset plot for APOE E3/E4 vs E3/E3
# Get colors
celltype_colors_path <- "path/to/colors.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Modify dataframe to duplicate colors for each diagnosis
celltype_colors <- rbind(celltype_colors, celltype_colors)
celltype_colors$predicted.celltype.l2[1:31] <- paste0("HC_", celltype_colors$predicted.celltype.l2[1:31])
celltype_colors$predicted.celltype.l2[32:62] <- paste0("AD_", celltype_colors$predicted.celltype.l2[32:62])
celltype_colors$predicted.celltype.l2 <- gsub(" ", "_",celltype_colors$predicted.celltype.l2)

# Initialize list
sig_genes_ls <- list()
for (diag in diag_list) {
  diag <- gsub(" ", "_", diag)
  print(diag)
  for (cluster in t_cell_clusters) {
    cluster_name <- gsub(" ", "_", cluster)
    print(cluster_name)
    # Load in degs
    tryCatch({
      degs <- read.csv(paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe34_vs_33_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      print(cluster)
      print(length(sig_genes[,1]))
      print(head(sig_genes))
      
      # Add sig genes to list
      sig_genes_ls[[paste0(diag, "_", cluster_name)]] <- sig_genes$X
    }, error = function(e) {
      NULL
    })
  }
}

sig_genes_ls

# Add groups of interest to dataframe for simple barplot (added to script 12-20-22)
# CD8 TCM
hc_cd8_tcm_34v33 <- length(sig_genes_ls[["Healthy_Control_CD8_TCM"]])
ad_cd8_tcm_34v33 <- length(sig_genes_ls[["Alzheimers_Disease_CD8_TCM"]])
df_CD8_TCM <- data.frame(hc_cd8_tcm_34v33,ad_cd8_tcm_34v33)
# CD4 CTL
hc_cd4_ctl_34v33 <- length(sig_genes_ls[["Healthy_Control_CD4_CTL"]])
ad_cd4_ctl_34v33 <- length(sig_genes_ls[["Alzheimers_Disease_CD4_CTL"]])
df_CD4_CTL <- data.frame(hc_cd4_ctl_34v33,ad_cd4_ctl_34v33)
# MAIT
hc_mait_34v33 <- length(sig_genes_ls[["Healthy_Control_MAIT"]])
ad_mait_34v33 <- length(sig_genes_ls[["Alzheimers_Disease_MAIT"]])
df_mait <- data.frame(hc_mait_34v33,ad_mait_34v33)

# Continue with upset
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]
names(sig_genes_ls) <- gsub("Healthy_Control", "HC", names(sig_genes_ls))
names(sig_genes_ls) <- gsub("Alzheimers_Disease", "AD", names(sig_genes_ls))

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "upset_apoe34_vs_33_sex_covar_genotype_comparisons_p0.01_lfc0.125.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()

# B.) UPset plot for APOE E4/E4 vs E3/E4
# Get colors
celltype_colors_path <- "path/to/colors.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Modify dataframe to duplicate colors for each diagnosis
celltype_colors <- rbind(celltype_colors, celltype_colors)
celltype_colors$predicted.celltype.l2[1:31] <- paste0("HC_", celltype_colors$predicted.celltype.l2[1:31])
celltype_colors$predicted.celltype.l2[32:62] <- paste0("AD_", celltype_colors$predicted.celltype.l2[32:62])
celltype_colors$predicted.celltype.l2 <- gsub(" ", "_",celltype_colors$predicted.celltype.l2)
# Initialize list
sig_genes_ls <- list()
for (diag in diag_list) {
  diag <- gsub(" ", "_", diag)
  print(diag)
  for (cluster in t_cell_clusters) {
    cluster_name <- gsub(" ", "_", cluster)
    print(cluster_name)
    # Load in degs
    tryCatch({
      degs <- read.csv(paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_34_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      print(cluster)
      print(length(sig_genes[,1]))
      print(head(sig_genes))
      
      # Add sig genes to list
      sig_genes_ls[[paste0(diag, "_", cluster_name)]] <- sig_genes$X
    }, error = function(e) {
      NULL
    })
  }
}

sig_genes_ls

# Add groups of interest to dataframe for simple barplot (added to script 12-20-22)
# CD8 TCM
hc_cd8_tcm_44v34 <- length(sig_genes_ls[["Healthy_Control_CD8_TCM"]])
ad_cd8_tcm_44v34 <- length(sig_genes_ls[["Alzheimers_Disease_CD8_TCM"]])
df_CD8_TCM <- cbind(df_CD8_TCM, data.frame(hc_cd8_tcm_44v34,ad_cd8_tcm_44v34))
# CD4 CTL
hc_cd4_ctl_44v34 <- length(sig_genes_ls[["Healthy_Control_CD4_CTL"]])
ad_cd4_ctl_44v34 <- length(sig_genes_ls[["Alzheimers_Disease_CD4_CTL"]])
df_CD4_CTL <- cbind(df_CD4_CTL, data.frame(hc_cd4_ctl_44v34,ad_cd4_ctl_44v34))
# MAIT
hc_mait_44v34 <- length(sig_genes_ls[["Healthy_Control_MAIT"]])
ad_mait_44v34 <- length(sig_genes_ls[["Alzheimers_Disease_MAIT"]])
df_mait <- cbind(df_mait, data.frame(hc_mait_44v34,ad_mait_44v34))

# Continue with upset
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]
names(sig_genes_ls) <- gsub("Healthy_Control", "HC", names(sig_genes_ls))
names(sig_genes_ls) <- gsub("Alzheimers_Disease", "AD", names(sig_genes_ls))

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "upset_apoe44_vs_34_sex_covar_genotype_comparisons_p0.01_lfc0.125.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()

# C.) UPset plot for APOE E4/E4 vs E3/E3
# Get colors
celltype_colors_path <- "path/to/colors.csv"
celltype_colors <- read.csv(celltype_colors_path)
celltype_colors[which(celltype_colors$predicted.celltype.l2 == "NK CD56bright"), "predicted.celltype.l2"] <- "NK_CD56bright"

# Modify dataframe to duplicate colors for each diagnosis
celltype_colors <- rbind(celltype_colors, celltype_colors)
celltype_colors$predicted.celltype.l2[1:31] <- paste0("HC_", celltype_colors$predicted.celltype.l2[1:31])
celltype_colors$predicted.celltype.l2[32:62] <- paste0("AD_", celltype_colors$predicted.celltype.l2[32:62])
celltype_colors$predicted.celltype.l2 <- gsub(" ", "_",celltype_colors$predicted.celltype.l2)
# Initialize list
sig_genes_ls <- list()
for (diag in diag_list) {
  diag <- gsub(" ", "_", diag)
  print(diag)
  for (cluster in t_cell_clusters) {
    cluster_name <- gsub(" ", "_", cluster)
    print(cluster_name)
    # Load in degs
    tryCatch({
      degs <- read.csv(paste0(output_dir, diag, "_", cluster_name, "_clonal_apoe44_vs_33_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      print(cluster)
      print(length(sig_genes[,1]))
      print(head(sig_genes))
      
      # Add sig genes to list
      sig_genes_ls[[paste0(diag, "_", cluster_name)]] <- sig_genes$X
    }, error = function(e) {
      NULL
    })
  }
}

sig_genes_ls

# Add groups of interest to dataframe for simple barplot (added to script 12-20-22)
# CD8 TCM
hc_cd8_tcm_44v33 <- length(sig_genes_ls[["Healthy_Control_CD8_TCM"]])
ad_cd8_tcm_44v33 <- length(sig_genes_ls[["Alzheimers_Disease_CD8_TCM"]])
df_CD8_TCM <- cbind(df_CD8_TCM, data.frame(hc_cd8_tcm_44v33,ad_cd8_tcm_44v33))
# CD4 CTL
hc_cd4_ctl_44v33 <- length(sig_genes_ls[["Healthy_Control_CD4_CTL"]])
ad_cd4_ctl_44v33 <- length(sig_genes_ls[["Alzheimers_Disease_CD4_CTL"]])
df_CD4_CTL <- cbind(df_CD4_CTL, data.frame(hc_cd4_ctl_44v33,ad_cd4_ctl_44v33))
# MAIT
hc_mait_44v33 <- length(sig_genes_ls[["Healthy_Control_MAIT"]])
ad_mait_44v33 <- length(sig_genes_ls[["Alzheimers_Disease_MAIT"]])
df_mait <- cbind(df_mait, data.frame(hc_mait_44v33,ad_mait_44v33))

# Continue with upset
sig_genes_ls <- sig_genes_ls[ lapply(sig_genes_ls, length) > 0 ]
names(sig_genes_ls) <- gsub("Healthy_Control", "HC", names(sig_genes_ls))
names(sig_genes_ls) <- gsub("Alzheimers_Disease", "AD", names(sig_genes_ls))

# Find number of genes in each set
num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
for (key in names(sig_genes_ls)) {
  num_degs <- rbind(num_degs, data.frame(key, length(sig_genes_ls[[key]])))
}

# Get order of sets for coloring
num_degs <- num_degs[order(-num_degs[, 2]),]

celltype_colors <- celltype_colors[match(num_degs[,1], celltype_colors$predicted.celltype.l2),]
celltype_colors <- celltype_colors[ which(celltype_colors$predicted.celltype.l2 %in% names(sig_genes_ls)),]

# Create upset plot and export
pdf(file = paste0(output_dir, "upset_apoe44_vs_33_sex_covar_genotype_comparisons_p0.01_lfc0.125.pdf"))
upset(
  fromList(sig_genes_ls),
  nsets = length(sig_genes_ls),
  order.by = "freq",
  sets.bar.color = celltype_colors$color
)
dev.off()

# Make simple bar plots for comparisons of interest

# CD8 TCM
df_CD8_TCM <- t(df_CD8_TCM)
df_CD8_TCM <- data.frame(cbind(rownames(df_CD8_TCM), df_CD8_TCM))
df_CD8_TCM$X1 <- factor(df_CD8_TCM$X1, levels = df_CD8_TCM$X1)
df_CD8_TCM$X2 <- as.numeric(df_CD8_TCM$X2)
df_CD8_TCM$X3 <- as.character(celltype_colors$color[13])
x <- ggplot(df_CD8_TCM, aes(x = X1, y = X2, fill = "#ffd700") ) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rep(celltype_colors$color[13], 6))
set_panel_size(x,
               file = paste0(output_dir, "clonal_CD8_TCM_num_degs.pdf"),
               width = unit(5, "in"), height = unit(4, "in")
)

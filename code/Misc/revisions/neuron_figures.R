# create figures for neuron revision

# BiocManager::install("Nebulosa")

suppressMessages({
  library(Signac)
  library(Seurat)
  library(tidyverse)
  library(Nebulosa)
})

rm(list = ls())
setwd("/path/to/directory/")
source("/path/to/helper_functions.R")

#### Genes Expressed over Accessible Regions ####
# identify object paths
path_RNA <- "/path/to/seurat/object"
path_ATAC_CD4 <- "/path/to/seurat/object"
path_ATAC_nonCD4 <- "/path/to/seurat/object"

# read in RNA data
load(path_RNA)

# average num features per cell type
gene_counts <- s@meta.data %>%
  as_tibble() %>%
  select(orig.ident, predicted.celltype.l2, nFeature_RNA) %>%
  group_by(orig.ident, predicted.celltype.l2) %>%
  summarize(mean_genes = mean(nFeature_RNA),
            sd_genes = sd(nFeature_RNA)) %>%
  dplyr::rename(cell_type = predicted.celltype.l2)

rm(s)
gc()

# read in ATAC CD4
s <- readRDS(path_ATAC_CD4)

# average num features per cell type
peak_counts_CD4 <- s@meta.data %>%
  as_tibble() %>%
  select(Sample, predictedGroup, nFeature_peaks) %>%
  group_by(Sample, predictedGroup) %>%
  summarize(mean_peaks = mean(nFeature_peaks),
            sd_peaks = sd(nFeature_peaks))

rm(s)
gc()

# read in ATAC non CD4
s <- readRDS(path_ATAC_nonCD4)

# average num features per cell type
peak_counts_nonCD4 <- s@meta.data %>%
  as_tibble() %>%
  select(Sample, predictedGroup, nFeature_peaks) %>%
  group_by(Sample, predictedGroup) %>%
  summarize(mean_peaks = mean(nFeature_peaks),
            sd_peaks = sd(nFeature_peaks))

rm(s)
gc()

# combine
features_by_cell_type <- rbind(peak_counts_CD4, peak_counts_nonCD4) %>%
  group_by(predictedGroup) %>%
  summarize(mean_peaks = mean(mean_peaks)) %>%
  left_join(summarize(group_by(gene_counts, cell_type), mean_genes = mean(mean_genes)), by = c("predictedGroup" = "cell_type")) %>%
  left_join(select(color_map, predicted.celltype.l2, new_color), by = c("predictedGroup" = "predicted.celltype.l2")) %>%
  arrange(predictedGroup)

atac_rna_cor <- cor.test(features_by_cell_type$mean_peaks, features_by_cell_type$mean_genes)

# plot
p <- ggplot(data = features_by_cell_type, aes(x = mean_peaks, y = mean_genes, color = predictedGroup)) +
  geom_point() +
  annotate("text", x = 7800, y = 4000, label = paste0("r = ", round(atac_rna_cor$estimate, 2), "\np = ", round(atac_rna_cor$p.value, 2))) +
  theme_Publication_blank() +
  scale_x_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 1000)) +
  scale_y_continuous(limits = c(0, 4200), breaks = seq(0, 4000, 1000)) +
  scale_color_manual(values = features_by_cell_type$new_color) +
  labs(x = "Mean Number of Accessible Regions", y = "Mean Number of Genes Expressed", color = "Cell Type")

pdf(file = "./counts_of_genes_by_regions.pdf", width = 9, height = 7)
p
dev.off()


#### DARs over Accessible Regions ####
celltype_colors_path <- "/path/to/broad_celltype_color_map.csv"
color_map <- read.csv(celltype_colors_path)

da_base_dir <- "/path/to/da_broad_celltypes/"
input_dir <- "/path/to/da/results/"

comparison <- "ADvsHC"
padj.thresh <- 0.05
lfc.thresh <- 0.125

# Define cell types
cell_types <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells", "Other_T_Cells", "Other")

LR_dir <- paste0("/path/to/LR/results")

intersected_df <- read.csv(paste0(input_dir, "full_intersection/", comparison, "_",
                                  "padj", as.character(padj.thresh),
                                  "_lfc", as.character(lfc.thresh),
                                  "_LR+DESeq2_DARs.csv"), row.names = 1)

grab_gene_num <- function(cell_type, de_dir) {
  tryCatch({
    # Read in DEGs
    if (cell_type == "CD8+_T_Cells") {
      celltype_files <- list.files(de_dir, pattern = "CD8", full.names = TRUE)
    } else if (cell_type == "CD4+_T_Cells") {
      celltype_files <- list.files(de_dir, pattern = "CD4", full.names = TRUE)
    } else {
      celltype_files <- list.files(de_dir, pattern = cell_type, full.names = TRUE)
    }
    degs <- read.csv(grep(".csv", celltype_files, value = TRUE))
    return(nrow(degs))
  }, error = function(e) {
    return(0)
  })
}

# Get deg to gene ratio
gene_nums <- sapply(cell_types, grab_gene_num, de_dir = LR_dir, USE.NAMES = TRUE) %>% as.data.frame
names(gene_nums) <- "num_genes"
deg_nums <- apply(intersected_df, MARGIN = 2, function(x) {length(x[!is.na(x) & x != ""])}) %>% as.data.frame
rownames(deg_nums) <- rownames(gene_nums)
names(deg_nums) <- "num_degs"
gene_nums <- merge(gene_nums, deg_nums, by = 0, all.y = TRUE)
gene_nums$ratio <- gene_nums$num_degs / gene_nums$num_genes
gene_nums <- gene_nums[order(gene_nums$num_degs),,drop=FALSE]
save(gene_nums, file = "region_nums")

load("region_nums")

gene_nums <- gene_nums %>%
  left_join(color_map, by = c("Row.names" = "predicted.celltype.l2")) %>%
  filter(num_degs != 0) %>%
  arrange(Row.names)

degs_genes_cor <- cor.test(gene_nums$num_degs, gene_nums$num_genes)

# plot
p <- ggplot(data = gene_nums, aes(x = num_genes, y = num_degs, color = Row.names)) +
  geom_point() +
  annotate("text", x = 29800, y = 1000, label = paste0("r = ", round(degs_genes_cor$estimate, 2), "\np = ", round(degs_genes_cor$p.value, 2))) +
  theme_Publication_blank() +
  scale_x_continuous(limits = c(20000, 30000), breaks = seq(20000, 30000, 1000)) +
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 100)) +
  scale_color_manual(values = gene_nums$new_color) +
  labs(x = "Number of Regions Tested", y = "Number of DARs", color = "Cell Type")

pdf(file = "./DARs_by_region.pdf", width = 7, height = 6)
p
dev.off()


#### DEGs over Genes Expressed ####
celltype_colors_path <- "/path/to/celltype_color_map.csv"
color_map <- read.csv(celltype_colors_path)
color_map$predicted_celltype_l2 <- sapply(color_map$predicted.celltype.l2,
                                          function(x) {gsub(" ", "_", x)})

# Get deg to gene ratio
mast_dir <- paste0(de_base_dir, "/path/to/directory")
gene_nums <- sapply(cell_types, grab_gene_num, de_dir = mast_dir, USE.NAMES = TRUE) %>% as.data.frame
names(gene_nums) <- "num_genes"
deg_nums <- apply(intersected_df, MARGIN = 2, function(x) {length(x[!is.na(x) & x != ""])}) %>% as.data.frame
names(deg_nums) <- "num_degs"
gene_nums <- merge(gene_nums, deg_nums, by = 0, all.y = TRUE)
gene_nums$ratio <- gene_nums$num_degs / gene_nums$num_genes
gene_nums <- gene_nums[order(gene_nums$ratio),,drop=FALSE]

save(gene_nums, file = "gene_nums")

load("gene_nums")
gene_nums <- gene_nums %>%
  left_join(select(color_map, predicted_celltype_l2, new_color), by = c("Row.names" = "predicted_celltype_l2")) %>%
  filter(num_degs != 0)

degs_genes_cor <- cor.test(gene_nums$num_degs, gene_nums$num_genes)

# plot
p <- ggplot(data = gene_nums, aes(x = num_genes, y = num_degs, color = Row.names)) +
  geom_point() +
  annotate("text", x = 9800, y = 160, label = paste0("r = ", round(degs_genes_cor$estimate, 2), "\np = ", round(degs_genes_cor$p.value, 2))) +
  theme_Publication_blank() +
  scale_x_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 10)) +
  scale_color_manual(values = gene_nums$new_color) +
  labs(x = "Number of Genes Tested", y = "Number of DEGs", color = "Cell Type")

pdf(file = "./DEGs_by_genes.pdf", width = 9, height = 8)
p
dev.off()

#### Nebulosa Density UMAPs ####
# RNA UMAPs
seurat_object <- "/path/to/seurat/object"
load(seurat_object)
color_pal <- "inferno"

s$HC <- s$Diagnosis == "Healthy Control"
s$AD <- s$Diagnosis == "Alzheimers Disease"
p <- plot_density(s, c("HC", "AD"), pal = color_pal)
pdf(file = "./RNA_UMAP_Diagnosis.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$male <- s$Sex == "male"
s$female <- s$Sex == "female"
p <- plot_density(s, c("male", "female"), pal = color_pal)
pdf(file = "./RNA_UMAP_Sex.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$E3E3 <- s$APOE_genotype == "E3/E3"
s$E3E4 <- s$APOE_genotype == "E3/E4"
s$E4E4 <- s$APOE_genotype == "E4/E4"
p <- plot_density(s, c("E3E3", "E3E4", "E4E4"), pal = color_pal)
pdf(file = "./RNA_UMAP_APOE.pdf", width = 15, height = 4)
p + plot_layout(ncol = 3)
dev.off()

rm(s)
gc()



# ATAC UMAPs
path_ATAC_CD4 <- "path/to/seurat/umap/object"
path_ATAC_nonCD4 <- "path/to/seurat/umap/object"

# create separated ATAC UMAPs
s <- readRDS(path_ATAC_CD4)

s$HC <- s$Diagnosis == "Healthy Control"
s$AD <- s$Diagnosis == "Alzheimers Disease"
p <- plot_density(s, c("HC", "AD"), pal = color_pal)
pdf(file = "./ATAC_CD4_UMAP_Diagnosis.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$male <- s$Sex == "male"
s$female <- s$Sex == "female"
p <- plot_density(s, c("male", "female"), pal = color_pal)
pdf(file = "./ATAC_CD4_UMAP_Sex.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$E3E3 <- s$APOE_Genotype == "E3/E3"
s$E3E4 <- s$APOE_Genotype == "E3/E4"
s$E4E4 <- s$APOE_Genotype == "E4/E4"
p <- plot_density(s, c("E3E3", "E3E4", "E4E4"), pal = color_pal)
pdf(file = "./ATAC_CD4_UMAP_APOE.pdf", width = 15, height = 4)
p + plot_layout(ncol = 3)
dev.off()

rm(s, p)
gc()


s <- readRDS(path_ATAC_nonCD4)
s <- DietSeurat(s, dimreducs = "umap")

s$HC <- s$Diagnosis == "Healthy Control"
s$AD <- s$Diagnosis == "Alzheimers Disease"
p <- plot_density(s, c("HC", "AD"), pal = color_pal)
pdf(file = "./ATAC_nonCD4_UMAP_Diagnosis.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$male <- s$Sex == "male"
s$female <- s$Sex == "female"
p <- plot_density(s, c("male", "female"), pal = color_pal)
pdf(file = "./ATAC_nonCD4_UMAP_Sex.pdf", width = 10, height = 4)
p + plot_layout(ncol = 2)
dev.off()

s$E3E3 <- s$APOE_Genotype == "E3/E3"
s$E3E4 <- s$APOE_Genotype == "E3/E4"
s$E4E4 <- s$APOE_Genotype == "E4/E4"
p <- plot_density(s, c("E3E3", "E3E4", "E4E4"), pal = color_pal)
pdf(file = "./ATAC_nonCD4_UMAP_APOE.pdf", width = 15, height = 4)
p + plot_layout(ncol = 3)
dev.off()
rm(s, p)
gc()

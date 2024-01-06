# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 04-21-2023
# Written by: Natalie Piehl
# Summary: Calculate cell number per group
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
  library("car")
  library("multcomp")
})

# Organize inputs
cd4_atac_seurat_object <- "/path/to/object"
noncd4_atac_seurat_object <- "/path/to/object"
rna_seurat_object <-"/path/to/object"
output_dir <- "/path/to/directory/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Look at RNA
#-------------------------------------------------------------------------------

# Load object
load(rna_seurat_object)

# Isolate cell types, Diagnosis, and APOE
mt <- s[[c("Diagnosis", "APOE_genotype", "orig.ident")]]

# Change to more legible
mt$Diagnosis <- mapvalues(mt$Diagnosis, from = c("Healthy Control", "Alzheimers Disease"),
                          to = c("HC", "AD"))

# Add Diagnosis + APOE column
mt$Diagnosis_APOE <- paste(mt$Diagnosis, mt$APOE_genotype)
mt$Diagnosis_APOE <- factor(mt$Diagnosis_APOE, levels = c("HC E3/E3", "AD E3/E3",
                                                          "HC E3/E4", "AD E3/E4",
                                                          "HC E4/E4", "AD E4/E4"))

# Get cell type nums
cell_count <- table(mt[,c("orig.ident")]) %>% as.data.frame
cell_count <- merge(cell_count, unique(mt), by.x = "Var1", by.y = "orig.ident", all.x = TRUE)

# Generate boxplot
p <- ggplot(cell_count, aes(x = Diagnosis_APOE, y = Freq, fill = Diagnosis_APOE)) +
  geom_boxplot(fatten = 1, width = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 2) +
  theme_Publication_blank() +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  labs(y = "Cell Count per Sample",
       x = NULL,
       title = "RNA") +
  theme(text = element_text(size = 25),
        legend.title=element_blank(),
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 8))
p

# Export
set_panel_size(p, file = paste0(output_dir, "RNA_cell_num_per_sample_by_group.pdf"))
write.csv(cell_count, file = paste0(output_dir, "RNA_cell_num_per_sample_by_group.csv"))

# Run ANOVA
res.aov <- aov(Freq ~ Diagnosis_APOE, data = cell_count)

# Generate ANOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "RNA_cell_num_per_sample_by_group_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- TukeyHSD(res.aov)
write.csv(postHocs$Diagnosis_APOE, file = paste0(output_dir, "RNA_cell_num_per_sample_by_group_tukey_pvals.csv"))

#-------------------------------------------------------------------------------
# Look at ATAC
#-------------------------------------------------------------------------------

# Read in ATAC
rm(s)
s <- readRDS(noncd4_atac_seurat_object)

# Isolate cell types, Diagnosis, and APOE
mt <- s[[c("Diagnosis", "APOE_genotype", "Sample")]]

# Change to more legible
mt$Diagnosis <- mapvalues(mt$Diagnosis, from = c("Healthy Control", "Alzheimers Disease"),
                          to = c("HC", "AD"))

# Add Diagnosis + APOE column
mt$Diagnosis_APOE <- paste(mt$Diagnosis, mt$APOE_genotype)
mt$Diagnosis_APOE <- factor(mt$Diagnosis_APOE, levels = c("HC E3/E3", "AD E3/E3",
                                                          "HC E3/E4", "AD E3/E4",
                                                          "HC E4/E4", "AD E4/E4"))

# Get cell type nums
cell_count_noncd4 <- table(mt[,c("Sample")]) %>% as.data.frame
cell_count_noncd4 <- merge(cell_count_noncd4, unique(mt), by.x = "Var1", by.y = "Sample", all.x = TRUE)

# Read in ATAC
rm(s)
s <- readRDS(cd4_atac_seurat_object)

# Isolate cell types, Diagnosis, and APOE
mt <- s[[c("Diagnosis", "APOE_genotype", "Sample")]]

# Change to more legible
mt$Diagnosis <- mapvalues(mt$Diagnosis, from = c("Healthy Control", "Alzheimers Disease"),
                          to = c("HC", "AD"))

# Add Diagnosis + APOE column
mt$Diagnosis_APOE <- paste(mt$Diagnosis, mt$APOE_genotype)
mt$Diagnosis_APOE <- factor(mt$Diagnosis_APOE, levels = c("HC E3/E3", "AD E3/E3",
                                                          "HC E3/E4", "AD E3/E4",
                                                          "HC E4/E4", "AD E4/E4"))

# Get cell type nums
cell_count_cd4 <- table(mt[,c("Sample")]) %>% as.data.frame
cell_count_cd4 <- merge(cell_count_cd4, unique(mt), by.x = "Var1", by.y = "Sample", all.x = TRUE)

# Add together
cell_count <- cell_count_noncd4
cell_count$Freq <- cell_count$Freq + cell_count_cd4$Freq

# Generate boxplot
p <- ggplot(cell_count, aes(x = Diagnosis_APOE, y = Freq, fill = Diagnosis_APOE)) +
  geom_boxplot(fatten = 1, width = 0.5, outlier.shape = NA) +
  geom_jitter(shape = 21, position = position_jitter(0.2), size = 2) +
  theme_Publication_blank() +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  labs(y = "Cell Count per Sample",
       x = NULL,
       title = "ATAC") +
  theme(text = element_text(size = 25),
        legend.title=element_blank(),
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 8))
p

# Export
set_panel_size(p, file = paste0(output_dir, "ATAC_cell_num_per_sample_by_group.pdf"))
write.csv(cell_count, file = paste0(output_dir, "ATAC_cell_num_per_sample_by_group.csv"))


# Run ANOVA
res.aov <- aov(Freq ~ Diagnosis_APOE, data = cell_count)

# Generate ANOVA table
anova_table <- Anova(res.aov)
write.csv(anova_table, file = paste0(output_dir, "ATAC_cell_num_per_sample_by_group_anova_table.csv"))

# Run post-hoc test (Tukey)
postHocs <- TukeyHSD(res.aov)
write.csv(postHocs$Diagnosis_APOE, file = paste0(output_dir, "ATAC_cell_num_per_sample_by_group_tukey_pvals.csv"))
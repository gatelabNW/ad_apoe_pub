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
  library("ggridges")
})

# Organize inputs
ranges_path <- "/path/to/ranges/object"
cicero_base_dir <- "/path/to/cicero/results"
input_base_dir <- "/path/to/RNA+ATAC/correlation/results/"
output_dir <- "/path/to/output_dir"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Load in results
#-------------------------------------------------------------------------------

# Define cell type comparisons
comparisons <- list.dirs(input_base_dir, full.names = FALSE)[-c(1)]

# Define diagnoses
diagnoses <- c("", "Healthy_Control", "Alzheimers_Disease")

# Initialize lists
R_list <- list()
p_list <- list()
pcc_thresh_list <- c()
num_sig_conn_list <- c()
num_glCRE_list <- c()
num_CRElinkedgenes_list <- c()
type_ls <- list()

# Iterate through comparisons
for (comparison in comparisons) {
  # Iterate through diagnoses
  for (diagnosis in diagnoses) {
    tryCatch({
      # Read in results
      df <- readRDS(paste0(input_base_dir, comparison, "/", diagnosis, "__", comparison,
                           '_peak_gene_correlation.rds'))
      
      # Append R and p values to list
      R_list[[paste(comparison, diagnosis)]] <- data.frame(R = df$pcc, 
                                                           comparison = rep(comparison, nrow(df)),
                                                           diagnosis = rep(diagnosis, nrow(df)))
      p_list[[paste(comparison, diagnosis)]] <- data.frame(p = df$BH, 
                                                           comparison = rep(comparison, nrow(df)),
                                                           diagnosis = rep(diagnosis, nrow(df)))
      
      # Subset for 95th percentile of PCC and adjusted p val < .01
      pcc_thresh <- quantile(abs(df$pcc), 0.95) %>% as.vector()
      padj_thresh <- .05
      df <- df[which(abs(df$pcc) >= pcc_thresh & df$BH <= padj_thresh),]
      pcc_thresh_list <- c(pcc_thresh_list, pcc_thresh)
      
      # Calculate proportion of each peak type
      peak_types <- table(df$Peak2_type)
      
      # Create formatted
      peak_types_formatted <- list(Distal = 0,
                                   Exonic = 0,
                                   Intronic = 0,
                                   Promoter = 0)
      for (type in names(peak_types_formatted)) {
        if (type %in% names(peak_types)) {
          peak_types_formatted[[type]] <- peak_types[[type]]
        }
      }
      
      # Add proportions to list
      type_ls[[paste(comparison, diagnosis)]] <- as.vector(unlist(peak_types_formatted))
      
      # Add nums to vectors
      num_sig_conn_list <- c(num_sig_conn_list, nrow(df))
      num_glCRE_list <- c(num_glCRE_list, length(unique(df$Peak2)))
      num_CRElinkedgenes_list <- c(num_CRElinkedgenes_list, length(unique(df$Peak1_nearestGene)))
    }, error = function(e) {
      pcc_thresh_list <<- c(pcc_thresh_list, 0)
      num_sig_conn_list <<- c(num_sig_conn_list, 0)
      num_glCRE_list <<- c(num_glCRE_list, 0)
      num_CRElinkedgenes_list <<- c(num_CRElinkedgenes_list, 0)
      type_ls[[paste(comparison, diagnosis)]] <<- c(0,0,0,0)
    })
  }
}

#-------------------------------------------------------------------------------
# Visualize
# - ridge plot of distribution of R and p
#-------------------------------------------------------------------------------

# Generate R df
df <- do.call(rbind, R_list)
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparisons)

# Ridge plot of R values
p <- ggplot(df, aes(x = abs(R), y = comparison, fill = diagnosis)) +
  geom_density_ridges(alpha = 0.4) +
  theme_Publication_blank() +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("white", "gray", "red")) +
  theme(text = element_text(size = 20))
p

# Export
set_panel_size(p, file = paste0(output_dir, "correlation_coefficient_ridge_plot.pdf"),
               width = unit(4, "in"), height = unit(10, "in"))

# Generate p df
df <- do.call(rbind, p_list)
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparisons)

# Ridge plot of p values
p <- ggplot(df, aes(x = -log10(p), y = comparison, fill = diagnosis)) +
  geom_density_ridges(alpha = 0.5) +
  theme_Publication_blank() +
  scale_x_log10() +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("white", "gray", "red")) +
  theme(text = element_text(size = 20))
p

# Export
set_panel_size(p, file = paste0(output_dir, "adjusted_pval_ridge_plot.pdf"),
               width = unit(4, "in"), height = unit(10, "in"))

# Box plot of p values
p <- ggplot(df, aes(y = -log10(p), x = comparison, fill = diagnosis)) +
  geom_boxplot(alpha = 0.5) +
  theme_Publication_blank() +
  scale_y_log10() +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("white", "gray", "red")) +
  theme(text = element_text(size = 20))
p

# Export
set_panel_size(p, file = paste0(output_dir, "adjusted_pval_box_plot.pdf"),
               width = unit(10, "in"), height = unit(4, "in"))

#-------------------------------------------------------------------------------
# Visualize
# - heatmap of num sig correlations and pcc threshold
# - heatmap of num gl-CREs vs num Cre-linked genes
#-------------------------------------------------------------------------------

# Define comparison levels
df <- data.frame(num_sig_conn = num_sig_conn_list,
                 comparison = rep(comparisons, each = 3),
                 diagnosis = rep(diagnoses, length(comparisons)))
ad_df <- df[which(df$diagnosis == "Alzheimers_Disease"), ]
ad_df <- ad_df[order(ad_df$num_sig_conn),]
comparison_levels <- unique(ad_df$comparison)

# Generate pcc df
df <- data.frame(pcc_thresh = pcc_thresh_list,
                     comparison = rep(comparisons, each = 3),
                     diagnosis = rep(diagnoses, length(comparisons)))
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparison_levels)

# Generate plot
p <- ggplot(df, aes(y = comparison, x = diagnosis, fill=pcc_thresh)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(pcc_thresh))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top", labels = c("all", "HC", "AD")) +
  ggtitle("R threshold (95th percentile)") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "pcc_thresh_heatmap.pdf"),
               width = unit(3, "in"), height = unit(8, "in"))

# Generate num_sig_conn df
df <- data.frame(num_sig_conn = num_sig_conn_list,
                 comparison = rep(comparisons, each = 3),
                 diagnosis = rep(diagnoses, length(comparisons)))
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparison_levels)

# Generate plot
p <- ggplot(df, aes(y = comparison, x = diagnosis, fill=num_sig_conn)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(num_sig_conn))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top", labels = c("all", "HC", "AD")) +
  ggtitle("# Significant Correlations") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_sig_conn_heatmap.pdf"),
               width = unit(3, "in"), height = unit(8, "in"))

# Generate glCRE df
df <- data.frame(num_glCRE = num_glCRE_list,
                 comparison = rep(comparisons, each = 3),
                 diagnosis = rep(diagnoses, length(comparisons)))
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparison_levels)

# Generate plot
p <- ggplot(df, aes(y = comparison, x = diagnosis, fill=num_glCRE)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(num_glCRE))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top", labels = c("all", "HC", "AD")) +
  ggtitle("# gene-linked CREs") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_glCRE_heatmap.pdf"),
               width = unit(3, "in"), height = unit(8, "in"))

# Generate Cre-linked-gene df
df <- data.frame(num_CRElinkedgenes = num_CRElinkedgenes_list,
                 comparison = rep(comparisons, each = 3),
                 diagnosis = rep(diagnoses, length(comparisons)))
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparison_levels)

# Generate plot
p <- ggplot(df, aes(y = comparison, x = diagnosis, fill=num_CRElinkedgenes)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(num_CRElinkedgenes))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top", labels = c("all", "HC", "AD")) +
  ggtitle("# cCRE-linked genes") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_CRElinkedgenes_heatmap.pdf"),
               width = unit(3, "in"), height = unit(8, "in"))

#-------------------------------------------------------------------------------
# Visualize
# - bar plot of peak types
#-------------------------------------------------------------------------------

# Convert to dataframe
res <- data.frame(type_ls) %>% t() %>% as.data.frame
colnames(res) <- c("Distal", "Exonic", "Intronic", "Promoter")

# Format
res <- res / rowSums(res)
res$comparison <- rep(comparisons, each = 3)
res$diagnosis <- rep(c("all", diagnoses[2:3]), length(comparisons))

# Convert to long format
res_long <- res %>%
  pivot_longer(cols = -c(comparison, diagnosis),
               names_to = "Peak_Type",
               values_to = 'Proportion')

# Generate plot
for (diagnosis in c("all", diagnoses[2:3])) {
  p <- ggplot(res_long[which(res_long$diagnosis == diagnosis),], aes(x = comparison, y = Proportion, fill = Peak_Type)) +
    geom_bar(stat = "identity", color = "black") +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = c("plum1", "violetred1", "skyblue", "palegreen")) +
    theme_Publication_blank() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "right",
          text = element_text(size = 24)) +
    ggtitle(diagnosis)
  p
  
  # Export plot
  set_panel_size(p, file = paste0(output_dir, diagnosis, "_gl-CRE_peak-type_composition.pdf"))
}

#-------------------------------------------------------------------------------
# Visualize
# - LFC of AD over HC sig correlations
#-------------------------------------------------------------------------------

# Generate num_sig_conn df
df <- data.frame(num_sig_conn = num_sig_conn_list,
                 comparison = rep(comparisons, each = 3),
                 diagnosis = rep(diagnoses, length(comparisons)))
df[which(df$diagnosis == ""), "diagnosis"] <- "all"
df$diagnosis <- factor(df$diagnosis, levels = c("all", diagnoses[2:3]))
df$comparison <- factor(df$comparison, levels = comparison_levels)

# Subset for AD and HC
df <- df[which(df$diagnosis %in% diagnoses[2:3]),]

# Convert to wide format
df_wide <- df %>%
  pivot_wider(names_from = "diagnosis",
               values_from = 'num_sig_conn')

# Add AD over HC num column
df_wide$NumSigCorr_ADoverHC_FC = df_wide$Alzheimers_Disease / df_wide$Healthy_Control
df_wide$label = rep("label", nrow(df_wide))
df_wide <- df_wide[is.finite(df_wide$NumSigCorr_ADoverHC_FC),]

# Generate plot
p <- ggplot(df_wide, aes(y = comparison, x = label, fill=log2(NumSigCorr_ADoverHC_FC))) +
  geom_tile() +
  geom_text(aes(label=scales::comma(NumSigCorr_ADoverHC_FC))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top", labels = c("all", "HC", "AD")) +
  ggtitle("LFC of # AD sig corr / # HC sig corr") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_sig_corr_lfc_heatmap.pdf"),
               width = unit(3, "in"), height = unit(8, "in"))

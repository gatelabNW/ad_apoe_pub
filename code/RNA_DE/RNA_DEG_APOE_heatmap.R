# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 11-29-2022
# Written by: Natalie Piehl
# Summary: Generate heatmaps comparing Diagnosis DE on different APOE genotypes
#
#-------------------------------------------------------------------------------
# Initialization

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("pheatmap")
})

# Organize inputs
de_base_dir <- "/path/to/DE/results/"
output_dir <- "/path/to/output/folder/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define thresholds
padj.thresh <- 0.01
lfc.thresh <- 0.25

#-------------------------------------------------------------------------------
# Get number of degs for each celltype (inter diagnosis)

# Define folders to compare
comparisons <- c("diagnosis_33", "diagnosis_34", "diagnosis_44")
deg_list <- list()

# Get celltype list
celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
                        function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()

for (comparison in comparisons) {
  # Define de dir
  # comparison <- comparisons[1]
  de_dir <- paste0(de_base_dir, comparison, "/out_NP_09-28-2022_covarSex/")
  
  # Make list of csvs
  csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
  
  # Initialize deg vector
  deg_vector <- c()
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, cell_type, "_AD_vs_HC_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, length(sig_genes$X))
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  # deg_list[[comparison]] <- unique(deg_vector)
  deg_list[[comparison]] <- deg_vector
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- c("E3/E3", "E3/E4", "E4/E4")

# Remove zero rows
data <- data[which(rowSums(data) != 0),]

# Add celltype column
data$celltype <- rownames(data)

# Convert to long format
data_long <- pivot_longer(
  data,
  cols = `E3/E3`:`E4/E4`, 
  names_to = "APOE",
  values_to = "Num_Degs"
)

# Generate plot
p <- ggplot(data_long, aes(x = celltype, y = APOE,
                          fill = Num_Degs)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red") +
  theme_Publication_blank() +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  theme(legend.position = "right",
        legend.direction="vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "pbmcsRNA_ADvsHC_numberOfDegs.pdf"),
               width = unit(7, "in"), height = unit(2, "in"))

# Generate with hclust
data_heatmap <- t(data[,1:3])
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         border_color = NA,
         filename = paste0(output_dir, "pbmcsRNA_ADvsHC_numberOfDegs_clust.pdf"),
         width = 5,
         height = 3
)

#-------------------------------------------------------------------------------
# Get number of degs for each celltype (inter APOE)

# Define folders to compare
comparisons <- c("apoe34vs33", "apoe44vs33", "apoe44vs34")
deg_list <- list()

# Get celltype list
celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
                        function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()

for (comparison in comparisons) {
  # Define de dir
  # comparison <- comparisons[1]
  de_dir <- paste0(de_base_dir, comparison, "_ad", "/out_NP_11-10-2022_covarSex/")

  # Make list of csvs
  csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)

  # Initialize deg vector
  deg_vector <- c()
  
  # Specify APOE label
  if (comparison == "apoe34vs33") {
    label <- "34_vs_33"
  } else if (comparison == "apoe44vs33") {
    label <- "44_vs_33"
  } else { label <- "44_vs_34"}

  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_ADonly_degs.csv"))

      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]

      # Add sig genes to list
      deg_vector <- c(deg_vector, length(sig_genes$X))
    }, error = function(e) {
      NULL
    })
  }

  # Add genes to list
  # deg_list[[comparison]] <- unique(deg_vector)
  deg_list[[comparison]] <- deg_vector
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- comparisons
ad_data <- data

for (comparison in comparisons) {
  # Define de dir
  # comparison <- comparisons[1]
  de_dir <- paste0(de_base_dir, comparison, "_hc", "/out_NP_11-10-2022_covarSex/")
  
  # Make list of csvs
  csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
  
  # Initialize deg vector
  deg_vector <- c()
  
  # Specify APOE label
  if (comparison == "apoe34vs33") {
    label <- "34_vs_33"
  } else if (comparison == "apoe44vs33") {
    label <- "44_vs_33"
  } else { label <- "44_vs_34"}
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_HConly_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, length(sig_genes$X))
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  # deg_list[[comparison]] <- unique(deg_vector)
  deg_list[[comparison]] <- deg_vector
}

# Convert to dataframe
data <- as.data.frame(deg_list)
rownames(data) <- celltype_list
colnames(data) <- comparisons
hc_data <- data

# Find difference
diff <- ad_data / hc_data

# Calculate LFC
diff <- log2(diff)

# Remove nonfinite rows
diff <- diff[is.finite(rowSums(diff)),]

# Add celltype column
diff$celltype <- rownames(diff)

# Convert to long format
diff_long <- pivot_longer(
  diff,
  cols = `apoe34vs33`:`apoe44vs34`,
  names_to = "APOE",
  values_to = "Log2Fold_Change_ADoverHC"
)

# Generate plot
p <- ggplot(diff_long, aes(x = celltype, y = APOE,
                           fill = Log2Fold_Change_ADoverHC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_Publication_blank() +
  scale_x_discrete(position = 'top', expand = c(0,0)) +
  scale_y_discrete(limits = rev, expand = c(0,0)) +
  theme(legend.position = "right",
        legend.direction="vertical",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "pbmcsRNA_APOEcomp_FC-OfDegsADvsHC.pdf"),
               width = unit(8, "in"), height = unit(2, "in"))

# Generate with hclust
data_heatmap <- t(diff[,1:3])
colfunc <- colorRampPalette(c("blue", "white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         breaks=seq(-5, 5, length.out=101),
         border_color = NA,
         filename = paste0(output_dir, "pbmcsRNA_APOEcomp_FC-OfDegsADvsHC_clust.pdf"),
         width = 5,
         height = 3
)

#-------------------------------------------------------------------------------
# Get number of degs for each celltype (inter APOE, just num)

# Define folders to compare
comparisons <- c("apoe34vs33", "apoe44vs33", "apoe44vs34")
deg_list <- list()

# Get celltype list
celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
                        function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()

for (comparison in comparisons) {
  # Define de dir
  # comparison <- comparisons[1]
  de_dir <- paste0(de_base_dir, comparison, "_ad", "/out_NP_11-10-2022_covarSex/")
  
  # Make list of csvs
  csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
  
  # Initialize deg vector
  deg_vector <- c()
  
  # Specify APOE label
  if (comparison == "apoe34vs33") {
    label <- "34_vs_33"
  } else if (comparison == "apoe44vs33") {
    label <- "44_vs_33"
  } else { label <- "44_vs_34"}
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_ADonly_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, length(sig_genes$X))
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
ad_deg_list <- deg_list

for (comparison in comparisons) {
  # Define de dir
  de_dir <- paste0(de_base_dir, comparison, "_hc", "/out_NP_11-10-2022_covarSex/")
  
  # Make list of csvs
  csvs <- list.files(de_dir, pattern = "*.csv", full.names = TRUE)
  
  # Initialize deg vector
  deg_vector <- c()
  
  # Specify APOE label
  if (comparison == "apoe34vs33") {
    label <- "34_vs_33"
  } else if (comparison == "apoe44vs33") {
    label <- "44_vs_33"
  } else { label <- "44_vs_34"}
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, cell_type, "_", label, "_HConly_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, length(sig_genes$X))
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
hc_deg_list <- deg_list

# Convert to dataframe
data_ad <- as.data.frame(ad_deg_list)
data_hc <- as.data.frame(hc_deg_list)
colnames(data_hc) <- paste0("HC_", comparisons)
colnames(data_ad) <- paste0("AD_", comparisons)
data <- cbind(data_hc, data_ad)
data <- data[,c(1,4,3,6,2,5)]
rownames(data) <- celltype_list

# Generate with hclust
data_heatmap <- t(data)
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         # breaks=seq(-5, 5, length.out=101),
         border_color = NA,
         filename = paste0(output_dir, "pbmcsRNA_APOEcomp_AD_HC_numDegs_clust.pdf"),
         width = 5,
         height = 3
)

#-------------------------------------------------------------------------------
# Get number of degs for each celltype (inter APOE, just num, clonal only)

# Define folders to compare
comparisons <- c("apoe34_vs_33", "apoe44_vs_33", "apoe44_vs_34")
deg_list <- list()

# Get celltype list
celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
                        function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()

for (comparison in comparisons) {
  # Define de dir
  de_dir <- "/path/to/tcr/DE/results/"

  # Initialize deg vector
  deg_vector <- c()

  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, "Alzheimers_Disease_", cell_type, "_clonal_", comparison, "_degs.csv"))

      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]

      # Add sig genes to list
      deg_vector <- c(deg_vector, nrow(sig_genes))
    }, error = function(e) {
      NULL
    })
  }

  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
ad_deg_list <- deg_list

cell_types <- c()
for (comparison in comparisons) {
  # Define de dir
  de_dir <- "/path/to/tcr/DE/results/"
  
  # Initialize deg vector
  deg_vector <- c()
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, "Healthy_Control_", cell_type, "_clonal_", comparison, "_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, nrow(sig_genes))
      
      cell_types <- c(cell_types, cell_type)
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
hc_deg_list <- deg_list

# Convert to dataframe
data_ad <- as.data.frame(ad_deg_list)
data_hc <- as.data.frame(hc_deg_list)
colnames(data_hc) <- paste0("HC_", comparisons)
colnames(data_ad) <- paste0("AD_", comparisons)
data <- cbind(data_hc, data_ad)
data <- data[,c(1,4,3,6,2,5)]
rownames(data) <- cell_types[1:9]
data <- data[rowSums(data) > 0,]
 
# Generate with hclust
data_heatmap <- t(data)
colfunc <- colorRampPalette(c("white", "red"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = TRUE,
         color = colfunc(100),
         # breaks=seq(-5, 5, length.out=101),
         border_color = NA,
         filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_numDegs_clust_RowClust.pdf"),
         width = 5,
         height = 3
)

#-------------------------------------------------------------------------------
# Get number of degs for each celltype (inter APOE, lfc, clonal only)

# Define folders to compare
comparisons <- c("apoe34_vs_33", "apoe44_vs_33", "apoe44_vs_34")
deg_list <- list()

# Get celltype list
celltype_list <- sapply(list.files(paste0(de_base_dir, "diagnosis_33/out_NP_09-28-2022_covarSex/"), pattern = "*.csv", full.names = TRUE),
                        function(x) {gsub("_AD_vs_HC_degs.csv", "", basename(x))}) %>% as.vector()

for (comparison in comparisons) {
  # Define de dir
  de_dir <- "/path/to/tcr/DE/results/"
  
  # Initialize deg vector
  deg_vector <- c()
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, "Alzheimers_Disease_", cell_type, "_clonal_", comparison, "_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, nrow(sig_genes))
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
ad_deg_list <- deg_list

cell_types <- c()
for (comparison in comparisons) {
  # Define de dir
  de_dir <- "/path/to/tcr/DE/results/"
  
  # Initialize deg vector
  deg_vector <- c()
  
  for (cell_type in celltype_list) {
    # Load in degs
    tryCatch({
      # Load in degs
      degs <- read.csv(paste0(de_dir, "Healthy_Control_", cell_type, "_clonal_", comparison, "_degs.csv"))
      
      # Identify sig genes
      sig_genes <- degs[which(degs$BH < padj.thresh & abs(degs$avg_log2FC) > lfc.thresh),]
      
      # Add sig genes to list
      deg_vector <- c(deg_vector, nrow(sig_genes))
      
      cell_types <- c(cell_types, cell_type)
    }, error = function(e) {
      NULL
    })
  }
  
  # Add genes to list
  deg_list[[comparison]] <- deg_vector
}
hc_deg_list <- deg_list

# # Convert to dataframe
# data_ad <- as.data.frame(ad_deg_list)
# data_hc <- as.data.frame(hc_deg_list)
# colnames(data_hc) <- paste0("HC_", comparisons)
# colnames(data_ad) <- paste0("AD_", comparisons)
# data <- cbind(data_hc, data_ad)
# data <- data[,c(1,4,3,6,2,5)]
# rownames(data) <- cell_types[1:9]
# data <- data[rowSums(data) > 0,]

# Convert to dataframe
data_ad <- as.data.frame(ad_deg_list)
data_hc <- as.data.frame(hc_deg_list)
colnames(data_hc) <- paste0("HC_", comparisons)
colnames(data_ad) <- paste0("AD_", comparisons)
rownames(data_hc) <- cell_types[1:9]
rownames(data_ad) <- cell_types[1:9]

# Find difference
diff <- data_ad / data_hc

# Calculate LFC
diff <- log2(diff)

# Remove nonfinite rows
diff <- diff[is.finite(rowSums(diff)),]

# Generate with hclust
data_heatmap <- t(data)
colfunc <- colorRampPalette(c("white", "red"))
# colfunc <- colorRampPalette(c("white", "red", "blue"))
pheatmap(as.matrix(data_heatmap),
         cluster_rows = FALSE,
         color = colfunc(100),
         # breaks=seq(-5, 5, length.out=101),
         border_color = NA,
         filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_numDegs_clust_HCtop.pdf"),
         # filename = paste0(output_dir, "pbmcsRNA_Clonal_APOEcomp_AD_HC_LFC_clust_RowClust.pdf"),
         width = 5,
         height = 3
)
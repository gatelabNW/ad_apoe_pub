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
})

# Organize inputs
ranges_path <- "/path/to/ranges/object"
input_base_dir <- "/path/to/cicero/results"
output_dir <- "/path/to/output_dir/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#-------------------------------------------------------------------------------
# Set up
#-------------------------------------------------------------------------------

# Read in ranges
ranges <- readRDS(ranges_path)

# Load in STAR gene data
ref_path <- "/path/to/ensembl/Homo_sapiens.GRCh38.109.chr.gtf.gz"
gene_anno <- rtracklayer::readGFF(ref_path)

# rename some columns to match plotting requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

#-------------------------------------------------------------------------------
# Make conn number plot
#-------------------------------------------------------------------------------

# Define cell types
cell_types <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells", "Other_T_Cells", "Other")

# Initialize conn list
conn_ls <- list()

# Iterate over celltypes
for (cell_type in cell_types[1:6]) {
  tryCatch({
    print(cell_type)
    # Load in conns
    hc_conns <- readRDS(paste0(input_base_dir, "Healthy_Control/", cell_type, "_conns.rds"))
    ad_conns <- readRDS(paste0(input_base_dir, "Alzheimers_Disease/", cell_type, "_conns.rds"))
    all_conns <- readRDS(paste0(input_base_dir, "all/", cell_type, "_conns.rds"))
    
    # remove all connections below threshold:
    ccan_cutoff <- 0.01
    hc_conns <- subset(hc_conns, coaccess >= ccan_cutoff)
    ad_conns <- subset(ad_conns, coaccess >= ccan_cutoff)
    all_conns <- subset(all_conns, coaccess >= ccan_cutoff)
    
    # create vector of number of conns
    res <- c(nrow(all_conns)/2, nrow(hc_conns)/2, nrow(ad_conns)/2)

    # append to list
    conn_ls[[cell_type]] <- res
  }, error = function(e) {
    message(e)
  })
}

# Bind results
conn_df <- do.call(rbind, conn_ls) %>% data.frame()

# Add names
names(conn_df) <- c('all', 'HC', 'AD')
conn_df$cell_type <- rownames(conn_df)

# Convert to long
conn_long <- pivot_longer(conn_df,
                          cols = 1:3,
                          names_to = "Diagnosis",
                          values_to = "Num_Conns")
conn_long$Diagnosis <- factor(conn_long$Diagnosis, levels = c("all", "HC", "AD"))

# Generate plot
p <- ggplot(conn_long, aes(y = cell_type, x = Diagnosis, fill=Num_Conns)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(Num_Conns))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top") +
  ggtitle("# Connections >= 0.01 coaccessibility") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_conns_0.01_allvsHCvsAD.pdf"),
               width = unit(3, "in"), height = unit(4, "in"))

#-------------------------------------------------------------------------------
# Make promoter conn number plot
#-------------------------------------------------------------------------------

# Define cell types
cell_types <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells", "Other_T_Cells", "Other")

# Initialize conn list
conn_ls <- list()

# Iterate over celltypes
for (cell_type in cell_types[1:6]) {
  tryCatch({
    print(cell_type)
    # Load in conns
    hc_conns <- readRDS(paste0(input_base_dir, "Healthy_Control/", cell_type, "_conns.rds"))
    ad_conns <- readRDS(paste0(input_base_dir, "Alzheimers_Disease/", cell_type, "_conns.rds"))
    all_conns <- readRDS(paste0(input_base_dir, "all/", cell_type, "_conns.rds"))
    
    # remove all connections below threshold:
    ccan_cutoff <- 0.01
    hc_conns <- subset(hc_conns, coaccess >= ccan_cutoff)
    ad_conns <- subset(ad_conns, coaccess >= ccan_cutoff)
    all_conns <- subset(all_conns, coaccess >= ccan_cutoff)
    
    # merge with annotations
    ad_conns <- merge(ad_conns, ranges, by.x = "Peak1", by.y = "sitename", all.x = TRUE)
    names(ad_conns)[4:5] <- c("Peak1_nearestGene", "Peak1_peakType")
    hc_conns <- merge(hc_conns, ranges, by.x = "Peak1", by.y = "sitename", all.x = TRUE)
    names(hc_conns)[4:5] <- c("Peak1_nearestGene", "Peak1_peakType")
    all_conns <- merge(all_conns, ranges, by.x = "Peak1", by.y = "sitename", all.x = TRUE)
    names(all_conns)[4:5] <- c("Peak1_nearestGene", "Peak1_peakType")
    
    # get only links that are connected to a gene promoter region
    ad_conns <- subset(ad_conns, Peak1_peakType == 'Promoter')
    hc_conns <- subset(hc_conns, Peak1_peakType == 'Promoter')
    all_conns <- subset(all_conns, Peak1_peakType == 'Promoter')
    
    # create vector of number of conns
    res <- c(nrow(all_conns), nrow(hc_conns), nrow(ad_conns))
    
    # append to list
    conn_ls[[cell_type]] <- res
  }, error = function(e) {
    message(e)
  })
}

# Bind results
conn_df <- do.call(rbind, conn_ls) %>% data.frame()

# Add names
names(conn_df) <- c('all', 'HC', 'AD')
conn_df$cell_type <- rownames(conn_df)

# Convert to long
conn_long <- pivot_longer(conn_df,
                          cols = 1:3,
                          names_to = "Diagnosis",
                          values_to = "Num_Conns")
conn_long$Diagnosis <- factor(conn_long$Diagnosis, levels = c("all", "HC", "AD"))

# Generate plot
p <- ggplot(conn_long, aes(y = cell_type, x = Diagnosis, fill=Num_Conns)) +
  geom_tile() +
  geom_text(aes(label=scales::comma(Num_Conns))) +
  scale_fill_fermenter(palette = "Blues", direction=1) +
  scale_x_discrete(position = "top") +
  ggtitle("# Promoter Connections >= 0.01 coaccessibility") +
  xlab('') + ylab('') +
  theme_Publication_blank() +
  theme(
    axis.ticks=element_blank(),
    legend.position = "right"
  )
p

# Export plot
set_panel_size(p, file = paste0(output_dir, "num_promoter_conns_0.01_allvsHCvsAD.pdf"),
               width = unit(3, "in"), height = unit(4, "in"))

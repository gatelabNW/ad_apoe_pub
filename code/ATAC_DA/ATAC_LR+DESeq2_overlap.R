# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
#
# Date: 05-02-2023
# Written by: Natalie Piehl
# Summary: Consolidate DAR lists
#
#-------------------------------------------------------------------------------
# Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("UpSetR")
  library("ComplexHeatmap")
})

# Organize inputs
da_base_dir <- "/path/to/DA/results/"
celltype_colors_path <- "/path/to/celltype_color_map"
output_base_dir <- "/path/to/output_dir/"
dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# Define comparison of interests
comparison <- "ADvsHC"

# Define thresholds
padj.thresh <- 0.05
lfc.thresh <- 0.125

# Define threshold specific output dir
output_dir <- paste0(output_base_dir, comparison, "/", "padj", as.character(padj.thresh),
                     "_lfc", as.character(lfc.thresh), "/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_base_dir, "full_intersection"), showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Prep
#-------------------------------------------------------------------------------

# Define cell types
cell_types <- c("B_Cells", "Monocytes", "Dendritic_Cells", "CD4+_T_Cells",
                "CD8+_T_Cells", "NK_Cells", "Other_T_Cells", "Other")

# Define comparison deg dirs
LR_dir <- paste0(da_base_dir, "main/out_NP_02-06-2023/", comparison, "/")
pseudobulk_dir <- paste0(da_base_dir, "DElegate/out_NP_05-05-2023/", comparison, "/")

# Initialize dar lists
LR_up_list <- list()
LR_down_list <- list()
pseudobulk_up_list <- list()
pseudobulk_down_list <- list()
down_list <- list()
up_list <- list()
combined_list <- list()

# Define dar function
grab_dars <- function(cell_type, de_dir, direction) {
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

    if (de_dir == pseudobulk_dir) {
      # Identify sig genes
      up_degs <- degs[which(degs$pvalue < padj.thresh & degs$log_fc > lfc.thresh),]
      down_degs <- degs[which(degs$pvalue < padj.thresh & degs$log_fc < -lfc.thresh),]
      
      # Save genes to list
      if (direction == "up") {
        return(up_degs$feature)
      } else {return(down_degs$feature)}
    } else {
      # Identify sig genes
      up_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC > lfc.thresh),]
      down_degs <- degs[which(degs$BH < padj.thresh & degs$avg_log2FC < -lfc.thresh),]
      
      # Save genes to list
      if (direction == "up") {
        return(up_degs$sitename)
      } else {return(down_degs$sitename)}
    }
  }, error = function(e) {
    return(c())
  })
}

# Initialize LR and pseudobulk overlap list
LR_pseudobulk_list <- list()

#-------------------------------------------------------------------------------
# Generate upsets per cell type
#-------------------------------------------------------------------------------

for (cell_type in cell_types) {
  # Grab LR DARs
  LR_up_list[[cell_type]] <- grab_dars(cell_type, LR_dir, "up")
  LR_down_list[[cell_type]] <- grab_dars(cell_type, LR_dir, "down")
  
  # Grab DESeq2 DARs
  pseudobulk_up_list[[cell_type]] <- grab_dars(cell_type, pseudobulk_dir, "up")
  pseudobulk_down_list[[cell_type]] <- grab_dars(cell_type, pseudobulk_dir, "down")

  # Generate up and down lists
  down_list[[cell_type]] <- list("LR" = LR_down_list[[cell_type]],
                                 "Pseudobulk" = pseudobulk_down_list[[cell_type]])
  up_list[[cell_type]] <- list("LR" = LR_up_list[[cell_type]],
                                 "Pseudobulk" = pseudobulk_up_list[[cell_type]])

  # Define color map
  color_map <- data.frame(method = c("LR", "Pseudobulk"),
                          color = c("cyan3", "plum2"))

  # Retain only nonzero groups
  down_list[[cell_type]] <- down_list[[cell_type]][ lapply(down_list[[cell_type]], length) > 0 ]
  up_list[[cell_type]] <- up_list[[cell_type]][ lapply(up_list[[cell_type]], length) > 0 ]

  # Up DARs --------------------------------------------------------------------
  if (length(up_list[[cell_type]]) > 1) {
    tryCatch({
      # Save results
      intersected_degs <- intersect(up_list[[cell_type]][["LR"]], up_list[[cell_type]][["Pseudobulk"]])
      write.csv(intersected_degs, paste0(output_dir, comparison, "_", cell_type, "_Up-LR+DESeq2_DARs.csv"))
    }, error = function(e) {
      message(e)
    })
  }

  # Down DARs ------------------------------------------------------------------
  if (length(down_list[[cell_type]]) > 1) {
    tryCatch({
      # Save results
      intersected_degs <- intersect(down_list[[cell_type]][["LR"]], down_list[[cell_type]][["Pseudobulk"]])
      write.csv(intersected_degs, paste0(output_dir, comparison, "_", cell_type, "_Down-LR+DESeq2_DARs.csv"))
    }, error = function(e) {
      message(e)
    })
  }

  # Combined DARs --------------------------------------------------------------
  # Check for discordant regions
  discordant_genes <- union(intersect(down_list[[cell_type]][["LR"]], up_list[[cell_type]][["LR"]]),
                            intersect(down_list[[cell_type]][["Pseudobulk"]], up_list[[cell_type]][["Pseudobulk"]]))
  if (length(discordant_genes) != 0) {
    print(discordant_genes)
    next
  }

  # Make combined DEG lists
  combined_list[[cell_type]] <- list("LR" = union(down_list[[cell_type]][["LR"]], up_list[[cell_type]][["LR"]]),
                                     "Pseudobulk" = union(down_list[[cell_type]][["Pseudobulk"]], up_list[[cell_type]][["Pseudobulk"]]))

  # Retain only nonzero groups
  combined_list[[cell_type]] <- combined_list[[cell_type]][ lapply(combined_list[[cell_type]], length) > 0 ]

  if (length(combined_list[[cell_type]]) > 1) {
    tryCatch({
      # Find number of genes in each set
      num_degs <- data.frame(matrix(ncol = 2, nrow = 0))
      for (key in names(combined_list[[cell_type]])) {
        num_degs <- rbind(num_degs, data.frame(key, length(combined_list[[cell_type]][[key]])))
      }

      # Get order of sets for coloring
      num_degs <- num_degs[order(-num_degs[, 2]),]
      color_map <- color_map[match(num_degs[,1], color_map$method),]

      # Create upset plot and export
      pdf(file = paste0(output_dir, comparison, "_", cell_type, "_DEGoverlap_upset.pdf"))
      print(upset(
        fromList(combined_list[[cell_type]]),
        text.scale = 3,
        nsets = length(combined_list[[cell_type]]),
        order.by = "freq",
        sets.bar.color = color_map$color
      ))
      dev.off()

      # Add intersection to list
      LR_pseudobulk_list[[cell_type]] <- intersect(combined_list[[cell_type]][["LR"]],
                                                     combined_list[[cell_type]][["Pseudobulk"]])
    }, error = function(e) {
      message(e)
    })
  }
}

# Export intersected list
max_degs <- 0
for (cell_type in cell_types) {
  degs <- LR_pseudobulk_list[[cell_type]]
  if (length(degs) > max_degs) {
    max_degs <- length(degs) 
  }
}
col_list <- list()
for (cell_type in cell_types) {
  degs <- LR_pseudobulk_list[[cell_type]]
  col_list[[cell_type]] <- c(degs, rep(NA, max_degs - length(degs)))
}
intersected_df <- do.call(cbind, col_list) %>% data.frame
intersected_df[is.na(intersected_df)] <- ""
write.csv(intersected_df, paste0(output_base_dir, "full_intersection/", comparison, "_",
                                 "padj", as.character(padj.thresh),
                                 "_lfc", as.character(lfc.thresh),
                                 "_LR+DESeq2_DARs.csv"))

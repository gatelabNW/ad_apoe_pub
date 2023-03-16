# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        AD APOE Project                             -----
# -----                                                                    -----
# -----                           Gate Lab                                 -----
# -----                     Northwestern University                        -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
# 
# Date: 06-30-2022
# Written by: Abhi Ramakrishnan
# Summary: Run ArchR (dev branch) pipeline on scATAC data
# Note: 180 G mem, 24 hours, 16 threads
#-------------------------------------------------------------------------------
# 1) Initialization
#-------------------------------------------------------------------------------

# Load in libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("ArchR")
  library("doMC")
  library("qpdf")
  library("readxl")
  library("rhdf5")
})


# Set random seed
set.seed(123)

# Set number of threads
addArchRThreads(threads = 26)

# Create project directory
proj_dir <- "/path/to/project/directory/"

ifelse(!dir.exists(proj_dir),
       dir.create(proj_dir, showWarnings = FALSE, recursive = TRUE), FALSE)
# Set working dir to project dir
setwd(proj_dir)

# Generate directory for manual outputs outside of project folder
output_dir <- "path/to/output/directory/"

ifelse(!dir.exists(output_dir),
       dir.create(output_dir, showWarnings = FALSE, recursive = TRUE), FALSE)

# -------------------------------------------------------------------------------
# 2) Generate ArchR arrow files and Project
# -------------------------------------------------------------------------------
# Create sample list
sample_list_all <- c("list of samples")

# Make list of fragment paths
fragments_dir <- "path/to/fragments/files/"
inputFiles <- c()
for (s in sample_list_all) {
  inputFiles <- c(inputFiles, s = paste0(fragments_dir, s, "/outs/fragments.tsv.gz"))
}
names(inputFiles) <- sample_list_all

# Specify genome
addArchRGenome("hg38")

# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles

# Add doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

# Generate ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = proj_dir,
  copyArrows = FALSE
)
# -------------------------------------------------------------------------------
# 3) Preliminary quality control
# -------------------------------------------------------------------------------
# Load project
proj <- loadArchRProject(proj_dir)
proj

# Remove samples which are not integrating with others on UMAP
# Skip and come back to this step after making initial UMAPs

# `%!in%` <- Negate(`%in%`)
# proj <- proj[which(proj$Sample %!in% c("list of samples"))]
# proj

# Subset project to contain only cells with tsse >= 10
idxPass <- which(proj$TSSEnrichment >= 10)
cellsPass <- proj$cellNames[idxPass]
proj <- proj[cellsPass, ]
min(proj$TSSEnrichment)
proj

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)
print("Project has been filtered to retain tsse > 10")
# -------------------------------------------------------------------------------
# 4) Add sample metadata
# -------------------------------------------------------------------------------
# Load project
proj <- loadArchRProject(proj_dir)

# Read in final samples metadata
scATAC_metadata <- read_excel("/path/to/metadata.xlsx", sheet=2)

# Clean up metadata and subset to retain columns of interest
meta <- scATAC_metadata[, c("matching_ID", "Age", "Diagnosis", "APOE_genotype", "Sex", "RNA_seq")]
colnames(meta) <- c("Sample", "Age", "Diagnosis", "APOE_Genotype", "Sex", "RNA")
# Get cell barcodes and sample IDs for each cell in project
cell_barcodes <- proj$cellNames
per_cell_sample_ID <- proj$Sample
# merge into df
per_cell_meta <- cbind(cell_barcodes, per_cell_sample_ID)
colnames(per_cell_meta) <- c("cell_barcodes", "Sample")
# merge cell each element of metadata by the sample ID
# Age
tmp <- merge(per_cell_meta, meta, all.x= TRUE, sort = FALSE)

# find order of samples in metadata and in project
meta_order <- unique(tmp$Sample)
proj_order <- unique(proj$Sample)
order_bind <- rbind(meta_order, proj_order) %>% t()
# check whether order of samples 
all.equal(meta_order, proj_order)

# Add correct metadata
proj$Age <- tmp$Age
proj$Diagnosis <- tmp$Diagnosis
proj$APOE_Genotype <- tmp$APOE_Genotype
proj$Sex <- tmp$Sex
proj$RNA <- tmp$RNA

# Check cell-level metadata is correct
df <- getCellColData(proj, select = c("Sample", "Diagnosis", "APOE_Genotype", "Age", "Sex", "RNA"))
tmp2 <- unique(df)
nrow(tmp2) # must equal 50
rownames(tmp2) <- NULL
all_equal(tmp2, data.frame(meta[which(meta$Sample %in% tmp$Sample),])) # must be TRUE

# Save project
proj <- saveArchRProject(ArchRProj = proj,
                         outputDirectory = proj_dir,
                         load = TRUE)

# -------------------------------------------------------------------------------
# 5) Visualizing QC metrics
# -------------------------------------------------------------------------------

# Load project
proj <- loadArchRProject(proj_dir)

# Plot fragment size distribtions per sample
p2 <- plotFragmentSizes(
  ArchRProj = proj,
  groupBy = "Sample",
  chromSizes = getChromSizes(proj),
  maxSize = 750,
  pal = NULL,
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotFragmentSizes")
)

# Plot tss enrichment per sample
p3 <- plotTSSEnrichment(ArchRProj = proj)

# Combine pdfs of QC metrics
tsse <- c()
doublet <- c()
fragment <- c()
qc_dir <- paste0(proj_dir, "/QualityControl/")
for (s in sample_list_all) {
  tsse <- c(tsse, paste0(qc_dir, s, "/", s, "-TSS_by_Unique_Frags.pdf" ))
  doublet <- c (doublet, paste0(qc_dir, s, "/", s, "-Doublet-Summary.pdf" ))
  fragment <- c (fragment, paste0(qc_dir, s, "/", s, "-Fragment_Size_Distribution.pdf" ))
}
pdf_combine(tsse, output = paste0(output_dir, "tsse_all.pdf"))
pdf_combine(doublet, output = paste0(output_dir, "doublet_all.pdf"))
pdf_combine(fragment, output = paste0(output_dir, "fragment_all.pdf"))

# Extract cellColData from project and make dataframe
cellColDataSample <- data.frame(proj@cellColData$Sample)
cellColDataNFrags <- data.frame(proj@cellColData$nFrags)
cellColDataDiag <- data.frame(proj@cellColData$Diagnosis)
cellColDataApoe <- data.frame(proj@cellColData$APOE_Genotype)
cellColDataTsse <- data.frame(proj@cellColData$TSSEnrichment)
# Group cells by diagnosis and apoe genotype into six sample groups
c <- cbind(cellColDataDiag, cellColDataApoe)
c$proj.cellColData.Diag.Apoe <- paste(c$proj.cellColData.Diagnosis,
                                      c$proj.cellColData.APOE_Genotype,
                                      sep = " ")
sample_qc <- cbind(cellColDataSample, cellColDataNFrags, cellColDataTsse, c)
write.csv(sample_qc, file= paste0(output_dir, "sample_qc.csv"))

# Define experimental groups and colors
diag_apoe <- c("Healthy Control E3/E3",
               "Healthy Control E3/E4",
               "Healthy Control E4/E4",
               "Alzheimers Disease E3/E3",
               "Alzheimers Disease E3/E4",
               "Alzheimers Disease E4/E4")
# Make shades of each color for AD groups
e3e3_col <- col2rgb("#63B8FF")
e3e3_shade <- rgb(red=.70*e3e3_col[1], green=.70*e3e3_col[2],blue=.70*e3e3_col[3], maxColorValue=255)
e3e4_col <- col2rgb("#C1FFC1")
e3e4_shade <- rgb(red=.70*e3e4_col[1], green=.70*e3e4_col[2],blue=.70*e3e4_col[3], maxColorValue=255)
e4e4_col <- col2rgb("#EEAEEE")
e4e4_shade <- rgb(red=.70*e4e4_col[1], green=.70*e4e4_col[2],blue=.70*e4e4_col[3], maxColorValue=255)
diag_apoe_colors <- c("Healthy Control E3/E3"="#63B8FF",
               "Healthy Control E3/E4"="#C1FFC1",
               "Healthy Control E4/E4"="#EEAEEE" ,
               "Alzheimers Disease E3/E3"=e3e3_shade,
               "Alzheimers Disease E3/E4"=e3e4_shade,
               "Alzheimers Disease E4/E4"=e4e4_shade)
# Order groups on plot
sample_qc$proj.cellColData.Diag.Apoe <- factor(sample_qc$proj.cellColData.Diag.Apoe,
                                                   levels = diag_apoe)
# Source custom plotting parameters
source(file="/path/to/helper_functions.R")
source(file="/projects/p31535/abhi/AD_APOE/code/plot/violin_plot.R")

# Violin plot of nFrags per group
a <- ggplot(sample_qc, aes(x=proj.cellColData.Diag.Apoe,
                           y=proj.cellColData.nFrags,
                           fill=proj.cellColData.Diag.Apoe)) +
  geom_violin() +
  scale_fill_manual(values=diag_apoe_colors,
                    name="Experimental Group",
                    guide = guide_legend(override.aes = list(shape= NA))) +
  labs(title = "Median Number of Fragments", x= NULL, y= "Number of Fragments/Cell") +
  custom_violin +
  theme_Publication_blank()
a
ggsave(a, file=paste0(output_dir, "nFrags_violin_group.pdf"))

# Violin plot of nFrags per group
a2 <- ggplot(sample_qc, aes(x=proj.cellColData.Sample,
                           y=proj.cellColData.nFrags,
                           fill=proj.cellColData.Sample)) +
  geom_violin() +
  scale_fill_manual(values=my_palette,
                    name= "Sample ID") +
  theme(legend.position = "none") +
  labs(title = "Median Number of Fragments", x= NULL, y= "Number of Fragments/Cell") +
  custom_violin
a2
ggsave(a2, dpi=300, height= 8, width = 4, units= "in", file=paste0(output_dir, "nFrags_violin_sample.pdf"))

# Violin plot of median tsse per group
b <- ggplot(sample_qc, aes(x=proj.cellColData.Diag.Apoe,
                           y=proj.cellColData.TSSEnrichment,
                           fill=proj.cellColData.Diag.Apoe)) +
  geom_violin() +
  scale_fill_manual(values=diag_apoe_colors,
                    name="Experimental Group",
                    guide = guide_legend(override.aes = list(shape= NA))) +
  labs(title = "Median TSS Enrichment", x= NULL, y= "TSS Enrichment Score") +
  custom_violin +
  theme_Publication_blank()
b
ggsave(b, file=paste0(output_dir, "tsse_violin_group.pdf"))
# Violin plot of median tsse per sample
b2 <- ggplot(sample_qc, aes(x=proj.cellColData.Sample,
                            y=proj.cellColData.TSSEnrichment,
                            fill=proj.cellColData.Sample)) +
  geom_violin() +
  scale_fill_manual(values=my_palette,
                    name= "Sample ID") +
  theme(legend.position = "none") +
  labs(title = "Median TSS Enrichment", x= NULL, y= "TSS Enrichment Score") +
  custom_violin
b2
ggsave(b2, dpi=300, height= 8, width = 4, units= "in", file=paste0(output_dir, "tsse_violin_sample.pdf"))

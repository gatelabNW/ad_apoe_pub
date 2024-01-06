# AD-APOE
# Abhirami Ramakrishnan
# #
# Biomarker boxplots for R4 Neuron revisions

# Load libraries
library(readxl)
library(tidyverse)
library(ggplot2)

# Set output dir
out_dir <- "/path/to/directory/"

# Read in CSF and plasma data
bio <- read.csv("/path/to/metadata.csv")

# Keep necessary columns
bio <- bio %>% 
  select(Diagnosis, PID_yr, CSF_pTau181_pg_mL,
         Plasma_pTau181_pg_mL, CSF_AB42_pg_mL,
         CSF_AB40_pg_mL, CSF_AB42_AB40_ratio)
  
# Remove outlier pTau entry
bio <- bio[!bio$PID_yr=="965_YR3", ]

# Create boxplots for each biomarker
# CSF_ptau
pdf(file = paste0(out_dir, "csf_ptau.pdf" )  )
ggplot(data=bio, mapping=aes(x=Diagnosis, y = CSF_pTau181_pg_mL, fill = Diagnosis)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  geom_jitter() +
  theme_classic()
dev.off()

# plasma_ptau
pdf(file = paste0(out_dir, "plasma_ptau.pdf" ))
ggplot(data=bio, mapping=aes(x=Diagnosis, y = Plasma_pTau181_pg_mL, fill = Diagnosis)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  geom_jitter() +
  theme_classic()
dev.off()
# csf_ab42
pdf(file = paste0(out_dir, "csf_ab42.pdf" ))
ggplot(data=bio, mapping=aes(x=Diagnosis, y = CSF_AB42_pg_mL, fill = Diagnosis)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  geom_jitter() +
  theme_classic()
dev.off()
# csf_ab40
pdf(file = paste0(out_dir, "csf_ab40.pdf" ))
ggplot(data=bio, mapping=aes(x=Diagnosis, y = CSF_AB40_pg_mL, fill = Diagnosis)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  geom_jitter() +
  theme_classic()
dev.off()
# csf_ab42_40_ratio
pdf(file = paste0(out_dir, "csf_ab42_40_ratio.pdf" ))
ggplot(data=bio, mapping=aes(x=Diagnosis, y = CSF_AB42_AB40_ratio, fill = Diagnosis)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  geom_jitter() +
  theme_classic()
dev.off()


wilcox.test(data = bio, CSF_AB42_AB40_ratio ~ Diagnosis )

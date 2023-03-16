#----------------------------------------------------------------
#----Calculate statistics for and plot age by sex, diagnosis, and APOE genotype
#----AD-APOE project
#----Abhi Ramakrishnan
#----Gate Lab
#----6.22.2022
#----------------------------------------------------------------
  
# Import libraries
library(tibble)
library(tidyverse) 
library(readxl) 
library(ggplot2)
library(RColorBrewer)
library(xlsx)
library(ggpubr)
library(rstatix)
library(car)

# Set working directory
setwd("path/to/directory/")

# Source helper functions
source("path/to/helper_functions.R")

# Define vector of sheet names
sheet_names <- c("scATAC_w_scRNA", "scATAC", "scRNA")
# Create empty dataframe for frequency counts
count_genotype_diagnosis <- data.frame(matrix(ncol=0, nrow=3))
rownames(count_genotype_diagnosis) <- c("APOE 3/3", "APOE 3/4", "APOE 4/4")

count_sex_diagnosis <- data.frame(matrix(ncol=0, nrow=2))
rownames(count_sex_diagnosis) <- c("Female", "Male")

count_sex_genotype_diagnosis <- data.frame(matrix(ncol=0, nrow=6))
rownames(count_sex_genotype_diagnosis) <- c("Female E3/E3",
                                   "Female E3/E4",
                                   "Female E4/E4",
                                   "Male E3/E3",
                                   "Male E3/E4",
                                   "Male E4/E4")

for (s in 1:3) {
  # Store sheet s in var t
  t <- read_excel("path/to/metadata.xlsx", sheet = s)
  # Remove unused patient samples
  if (s==2) {
    t <- t[!(t$Patient_ID %in% c("1111","696")),]
  } else if (s==3){
    t <- t[!(t$Patient_ID %in% c("1086","659")),]
  } else {
    t <- t[!(t$Patient_ID %in% c("1086","659","1111","696")),]
  }
  print(t$Patient_ID)
  # Set order of levels in Diagnosis variable to display HC first on plots
  t$Diagnosis <- factor(t$Diagnosis, levels = c("Healthy Control",
                                                        "Alzheimers Disease"))
  # Shapiro-Wilk test for normality
  a <- shapiro_test(t[which(t$Diagnosis == "Healthy Control" & t$APOE_genotype == "E3/E3"),], vars = "Age")
  b <- shapiro_test(t[which(t$Diagnosis == "Healthy Control" & t$APOE_genotype == "E3/E4"),], vars = "Age")
  c <- shapiro_test(t[which(t$Diagnosis == "Healthy Control" & t$APOE_genotype == "E4/E4"),], vars = "Age")
  d <- shapiro_test(t[which(t$Diagnosis == "Alzheimers Disease" & t$APOE_genotype == "E3/E3"),], vars = "Age")
  e <- shapiro_test(t[which(t$Diagnosis == "Alzheimers Disease" & t$APOE_genotype == "E3/E4"),], vars = "Age")
  f <- shapiro_test(t[which(t$Diagnosis == "Alzheimers Disease" & t$APOE_genotype == "E4/E4"),], vars = "Age")
  # combine results into table
  shapiro <- rbind(a[,2:3],b[,2:3],c[,2:3],d[,2:3],e[,2:3],f[,2:3])
  colnames(shapiro) <- c("Shapiro statistic", "p-val")
  rownames(shapiro) <- c("HC 3/3", "HC 3/4", "HC 4/4", "AD 3/3", "AD 3/4", "AD 4/4")
  # Anova multiple comparisons test for Age
  # Anova test; check F values to see if there is significant difference in variance of age by diagnosis or genotype
  anova <- anova_test(t, Age ~ Diagnosis * APOE_genotype)
  # Levene test to test for significant difference in variance between groups
  # CHANGE THE FOLLOWING LINE--transformation of age variable
  levene <- leveneTest(sqrt(Age) ~ Diagnosis*APOE_genotype, t)
  # Anova multiple comparisons test for Age
  tukey <- (t %>% tukey_hsd(Age ~ Diagnosis*APOE_genotype))
  # Create stats folder
  ifelse(!dir.exists( "outs/stats/"),
         dir.create("stats/", recursive = TRUE),
         print("Directory exists"))
  # Write stats results into csv
  write.xlsx(shapiro, file=file.path("outs/stats/", paste0(sheet_names[s],"-age_demographic_stats_filtered.xlsx")),
             sheetName="Shapiro-Wilk Test", row.names=TRUE)
  write.xlsx(anova, file=file.path("outs/stats/", paste0(sheet_names[s],"-age_demographic_stats_filtered.xlsx")),
             sheetName="2-way Anova", append = TRUE, row.names=TRUE)
  write.xlsx(levene, file=file.path("outs/stats/", paste0(sheet_names[s],"-age_demographic_stats_filtered.xlsx")),
             sheetName="Levene's test", append = TRUE, row.names=TRUE)
  write.xlsx(tukey, file=file.path("outs/stats/", paste0(sheet_names[s],"-age_demographic_stats_filtered.xlsx")),
             sheetName="Tukey Multiple Comparisons", append = TRUE, row.names=TRUE)
  # Plot age_v_diagnosis
  age_v_diagnosis <- ggplot(t, aes(fill=APOE_genotype, y=Age, x=Diagnosis)) +
    geom_boxplot() +
    theme_Publication_blank() +
    ggtitle(paste0("Age by Diagnosis for ", sheet_names[[s]])) +
    scale_fill_brewer(palette = "YlOrRd") +
    geom_dotplot(binaxis = "y",
                 stackdir = "center",
                 dotsize = 0.5,
                 position = position_dodge(0.75))
  show(age_v_diagnosis)
  ggsave(age_v_diagnosis, filename = paste0("outs/demographic_plots/", sheet_names[[s]], "_age_v_diagnosis_filtered.pdf"))
  # Plot age_v_genotype
  age_v_genotype <- ggplot(t, aes(fill=Diagnosis, y=Age, x=APOE_genotype)) +
    geom_boxplot() +
    theme_Publication_blank() +
    ggtitle(paste0("Age by Genotype and Diagnosis for ", sheet_names[[s]])) +
    scale_fill_brewer(palette = "PuRd") +
    geom_dotplot(binaxis = "y",
                 stackdir = "center",
                 dotsize = 0.5,
                 position = position_dodge(0.75))
  show(age_v_genotype)
  ggsave(age_v_genotype, filename = paste0("outs/demographic_plots/", sheet_names[[s]], "_age_v_genotype_filtered.pdf"))
  # Plot age_v_sex_diagnosis
  age_v_sex_diagnosis <- ggplot(t, aes(fill=Sex, y=Age, x=Diagnosis)) +
    geom_boxplot() +
    theme_Publication_blank() +
    ggtitle(paste0("Age by Sex and Diagnosis for ", sheet_names[[s]])) +
    scale_fill_brewer(c("ccebc5","7bccc4")) +
    geom_dotplot(binaxis = "y",
                 stackdir = "center",
                 dotsize = 0.5,
                 position = position_dodge(0.75))
  show(age_v_sex_diagnosis)
  ggsave(age_v_sex_diagnosis, filename = paste0("outs/demographic_plots/", sheet_names[[s]], "_age_v_sex_diagnosis_filtered.pdf"))
  # Count frequencies for each group:
  # Genotype and Diagnosis
  h <- t[which(t$Diagnosis == "Healthy Control"),] %>%
    count(APOE_genotype)
  d <- t[which(t$Diagnosis == "Alzheimers Disease"),] %>%
    count(APOE_genotype)
  count_genotype_diagnosis <- cbind(count_genotype_diagnosis, h[,2], d[,2])
  # Sex and Diagnosis
  h <- t[which(t$Diagnosis == "Healthy Control"),] %>%
    count(Sex)
  d <- t[which(t$Diagnosis == "Alzheimers Disease"),] %>%
    count(Sex)
  count_sex_diagnosis <- cbind(count_sex_diagnosis, h[,2], d[,2])
  # Genotype, Sex, and Diagnosis
  h <- t[which(t$Diagnosis == "Healthy Control"),] %>%
    count(Sex, APOE_genotype)
  d <- t[which(t$Diagnosis == "Alzheimers Disease"),] %>%
    count(Sex, APOE_genotype)
  count_sex_genotype_diagnosis <- cbind(count_sex_genotype_diagnosis, h[,3], d[,3])

}
# Add column names to frequency tables
colnames(count_genotype_diagnosis) <- c("common_HC", "common_AD", "scATAC_HC", "scATAC_AD", "scRNA_HC", "scRNA_AD")
colnames(count_sex_diagnosis) <- c("common_HC", "common_AD", "scATAC_HC", "scATAC_AD", "scRNA_HC", "scRNA_AD")
colnames(count_sex_genotype_diagnosis) <- c("common_HC", "common_AD", "scATAC_HC", "scATAC_AD", "scRNA_HC", "scRNA_AD")

# Export all frequency tables into an excel spreadhseet
write.xlsx(count_genotype_diagnosis, file="outs/demographic_frequencies_filtered.xlsx",
           sheetName="genotype_diagnosis", row.names=TRUE)
write.xlsx(count_sex_diagnosis, file="outs/demographic_frequencies_filtered.xlsx",
           sheetName="sex_diagnosis", append = TRUE, row.names=TRUE)
write.xlsx(count_sex_genotype_diagnosis, file="outs/demographic_frequencies_filtered.xlsx",
           sheetName="outs/sex_genotype_diagnosis", append = TRUE, row.names=TRUE)



# Import data for Run2 of the PredictTB samples
# 9/18/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction

# DuffyTools
library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


# cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
# cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
# cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbPalette2 <-  c("#bfbfbf", "#56B4E9")
# cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442","#CC79A7")
# c25 <- c(
#   "dodgerblue2", "#E31A1C", "green4",
#   "#6A3D9A","#FF7F00","black", "gold1",
#   "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
#   "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
# )
# c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4")
# c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4")


# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############

## PREDICTTB_RUN1
# This has been edited to include more metadata!
Run1_pipeSummary <- read.csv("Data/PredictTB_Run1/Pipeline.Summary.Details.csv") 

# Just get the samples I want 
Run1_pipeSummary <- Run1_pipeSummary %>% 
  filter(Type %in% c("THP1 spiked", "Week 0 sputum", "Week 2 sputum")) %>% 
  mutate(Run = "PredictTB_Run1")

# BROTH DATA FROM PROBETEST 5
ProbeTest5_pipeSummary <- read.csv("Data/ProbeTest5/ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") 
Broth_pipeSummary <- ProbeTest5_pipeSummary %>% 
  filter(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9", "THP1_1e6_1a_S28"))

## PREDICTTB_RUN2
# This has been edited to include more metadata!
Run2_pipeSummary <- read.csv("Data/PredictTB_Run2/Pipeline.Summary.Details.csv") 
Run2_pipeSummary <- Run2_pipeSummary %>%
  mutate(Run = "PredictTB_Run2")


# Merge the pipeSummaries
All_pipeSummary <- merge(Run1_pipeSummary, Run2_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, Broth_pipeSummary, all = T)

# Merge two columns
All_pipeSummary <- All_pipeSummary %>% mutate(Type = coalesce(Type, Sample_Type)) %>%
  mutate(Type2 = coalesce(Type2, Sample_Type)) %>% 
  mutate(Type2 = str_replace(Type2, "^THP1$", "THP1 spiked")) %>% 
  select(-Sample_Type) %>%
  select(-c(N_Splice, P_Splice, Replicate, RT, CFU_per_g, CFU_per_mL, Ra_cells, EukrRNADep, Hyb_Time, Probe, Probe_ng, Pooled_Set, X, mRNA_ng, ct, ttd)) %>%
  filter(SampleID != "Undetermined_S0")

# Make a second SampleID column
# Remove _S* From names
All_pipeSummary$SampleID2 <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)
All_pipeSummary$Run2 <- gsub(x = All_pipeSummary$Run, pattern = "PredictTB_", replacement = "")
All_pipeSummary <- All_pipeSummary %>% mutate(SampleID2 = paste0(Run2, "_", SampleID2))


###########################################################
########### LIST OF SAMPLES PASSING INSPECTION ############

GoodSampleList <- All_pipeSummary %>%
  filter(N_Genomic >= 1000000 & AtLeast.10.Reads >= (4499*0.75)) %>% # CHANGED TO BE 75% !!!!
  pull(SampleID2)

# Checked and all the names are unique

SputumSampleList <- GoodSampleList[grep("^W", GoodSampleList)]

BrothSampleList <- All_pipeSummary %>% 
  filter(str_detect(SampleID, "Broth")) %>%
  pull(SampleID)

GoodSamples_pipeSummary <- All_pipeSummary %>% filter(SampleID2 %in% GoodSampleList)

###########################################################
############### IMPORT AND PROCESS TPM VALUES #############

Run1_tpm <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.TPM.csv")
Run1_tpm <- Run1_tpm %>% select(X, THP1_1e6_1_S67, contains("W"))

# Just pull the tpm of the THP1 spiked from another run: THP1 1e6_1 (Predict rack 2 box 1 I04)
# Need THP1 1e6_1a from the Januaray run. Also need 
ProbeTest5_tpm <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") 
ProbeTest5_tpm_Broth <- ProbeTest5_tpm %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28)

Run2_tpm <- read.csv("Data/PredictTB_Run2/Mtb.Expression.Gene.Data.TPM.csv")
Run2_tpm <- Run2_tpm %>% select(-contains("Undetermined"))

# Adjust the names so they are slightly shorter
# names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# rownames(All_tpm) <- All_tpm[,1] # add the rownames

# Need to make sure there is a Gene column (gets lots)
# All_tpm <- All_tpm %>% 
#   rename(Gene = X) 


###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% 
  select(X, THP1_1e6_1_S67, contains("W")) %>%
  rename_with(~ paste0("Run1_", .), -X) # Add Run1 to the beginning of every column because some have the same name

ProbeTest5_RawReads <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.readsM_moreTrim.csv")
ProbeTest5_RawReads_Broth <- ProbeTest5_RawReads %>% 
  select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9, THP1_1e6_1a_S28) %>%
  rename_with(~ paste0("ProbeTest5_", .), -X) 

Run2_RawReads <- read.csv("Data/PredictTB_Run2/Mtb.Expression.Gene.Data.readsM.csv")
Run2_RawReads <- Run2_RawReads %>% 
  select(-contains("Undetermined")) %>% 
  rename_with(~ paste0("Run2_", .), -X) 
  

# Merge the RawReads I collected above
All_RawReads <- merge(Run1_RawReads, Run2_RawReads, all = T)
All_RawReads <- merge(All_RawReads, ProbeTest5_RawReads_Broth)

# Remove the _S at the end
names(All_RawReads) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_RawReads))

# # Just keep the samples passing filter
# All_RawReads <- All_RawReads %>% select("X", all_of(GoodSampleList), "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")
# 
# 
# 
###########################################################
################ MAKE TPM FROM Rv RAW READS ###############

# Bob's TPM includes all the non-coding RNAs, make a new TPM from just the protein-coding genes

# Keep only the protein coding Rv genes
All_RawReads_f <- All_RawReads %>%
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", X))

source("Function_CalculateTPM.R")
All_tpm_f <- CalculateTPM_RvOnly(All_RawReads_f)


# Remove the _S at the end
names(All_tpm_f) = gsub(pattern = "_S[0-9]+$", replacement = "", x = names(All_tpm_f))

# Subset just the GoodSamples
GoodSamples_tpmf <- All_tpm_f %>% select(all_of(GoodSampleList))



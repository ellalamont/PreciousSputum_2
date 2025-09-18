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
# cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442","#CC79A7")
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
  filter(SampleID %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9"))

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
  select(-Sample_Type) %>%
  select(-c(N_Splice, P_Splice, Replicate, RT, CFU_per_g, CFU_per_mL, Ra_cells, EukrRNADep, Hyb_Time, Probe, Probe_ng, Pooled_Set, X)) %>%
  filter(SampleID != "Undetermined_S0")


# Reorder things
# NOT DONE YET!
# All_pipeSummary$Drug <- as.character(All_pipeSummary$Drug)
# ordered_Drug <- c("Untreated", "RIF")
# All_pipeSummary$Drug <- factor(All_pipeSummary$Drug, levels = ordered_Drug)
# All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)


###########################################################
########### LIST OF SAMPLES PASSING INSPECTION ############

GoodSampleList <- All_pipeSummary %>%
  filter(N_Genomic >= 1000000 & AtLeast.10.Reads >= (4499*0.8)) %>% 
  pull(SampleID)

SputumSampleList <- GoodSampleList[grep("^W", GoodSampleList)]

BrothSampleList <- All_pipeSummary %>% 
  filter(str_detect(SampleID, "Broth")) %>%
  pull(SampleID)



# 9/18/25 BELOW IS NOT DONE YET!!!!

###########################################################
############### IMPORT AND PROCESS TPM VALUES #############

Run1_tpm <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.TPM.csv")
Run1_tpm <- Run1_tpm %>% select(-contains("TBAIT"))

# Just pull the tpm of the THP1 spiked from another run: THP1 1e6_1 (Predict rack 2 box 1 I04)
# Need THP1 1e6_1a from the Januaray run. Also need 
ProbeTest5_tpm <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") 
ProbeTest5_tpm_subset <- ProbeTest5_tpm %>% select(X, THP1_1e6_1a_S28)
ProbeTest5_tpm_Broth <- ProbeTest5_tpm %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9)

# Get the Marmoset TPM
ProbeTest3_tpm <- read.csv("Data/ProbeTest3/ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")
ProbeTest3_tpm_marm <- ProbeTest3_tpm %>% select(X, BQ12_10_Probe_3A_S29, BQ12_3_Probe_4A_50_S27, BQ12_8_Probe_4A_50_S28)




# Adjust the names so they are slightly shorter
# names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

# rownames(All_tpm) <- All_tpm[,1] # add the rownames

# Need to make sure there is a Gene column (gets lots)
# All_tpm <- All_tpm %>% 
#   rename(Gene = X) 


###########################################################
###### MAKE TPM WITH ALL CLINICAL AND ANIMAL MODELS #######

# Merge the tpms I collected above
All_tpm <- merge(Run1_tpm, ProbeTest3_tpm_marm, all = T)
All_tpm <- merge(All_tpm, ProbeTest5_tpm_Broth)

# Just keep the samples passing filter
All_tpm <- All_tpm %>% select("X", all_of(GoodSampleList), "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")



###########################################################
############### IMPORT AND PROCESS RAW READS ##############

Run1_RawReads <- read.csv("Data/PredictTB_Run1/Mtb.Expression.Gene.Data.readsM.csv")
Run1_RawReads <- Run1_RawReads %>% select(-contains("TBAIT"))

ProbeTest5_RawReads <- read.csv("Data/ProbeTest5/ProbeTest5_Mtb.Expression.Gene.Data.readsM_moreTrim.csv") 
ProbeTest5_RawReads_Broth <- ProbeTest5_RawReads %>% select(X, H37Ra_Broth_4_S7, H37Ra_Broth_5_S8, H37Ra_Broth_6_S9)

ProbeTest3_RawReads <- read.csv("Data/ProbeTest3/Mtb.Expression.Gene.Data.readsM_old.csv")
ProbeTest3_RawReads_marm <- ProbeTest3_RawReads %>% select(X, BQ12_10_Probe_3A_S29, BQ12_3_Probe_4A_50_S27, BQ12_8_Probe_4A_50_S28)

# Merge the RawReads I collected above
All_RawReads <- merge(Run1_RawReads, ProbeTest3_RawReads_marm, all = T)
All_RawReads <- merge(All_RawReads, ProbeTest5_RawReads_Broth)

# Just keep the samples passing filter
All_RawReads <- All_RawReads %>% select("X", all_of(GoodSampleList), "H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9")



###########################################################
################ NEW DATA WITH FEWER SPUTUM ###############

# Don't want to include all the sputum in the paper, want a random subset of ~12

# First remove all the relapse, then randomly choose
sputumOnly_pipeSummary <- All_pipeSummary %>% filter(str_detect(Type, "sputum"))
W0CureSputumOnly_pipeSummary <- sputumOnly_pipeSummary %>% filter(Outcome == "Cure" & Week == "Week 0")

# Randomly choose 12 sputum samples
set.seed(5) # 3, 9
SputumSubset_pipeSummary <- slice_sample(sputumOnly_pipeSummary, n = 12)
SputumSubset_list <- SputumSubset_pipeSummary %>% pull(SampleID)










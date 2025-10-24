# iModulons gene sets comparisons
# E. Lamont
# 10/24/25 Run1-3


# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result
# I know the UP and DOWN files are a little different, just working with the UP files for now, should maybe ask Bob if that is okay at some point

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=7),
        legend.title = element_text(size = 7),
        plot.title = element_text(size=7), 
        axis.title.x = element_text(size=7), 
        axis.text.x = element_text(angle = 0, size=7, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
        axis.title.y = element_text(size=7),
        axis.text.y = element_text(size=7), 
        plot.subtitle = element_text(size=7)# , 
        # plot.margin = margin(10, 10, 10, 20)# ,
  )

###########################################################
############ FUNCTION: ONE iMODULON CATEGORY ##############

BubblePlot_Function <- function(df, Category) {
  
  ## Generate a single bubble plot given a specific iModulon category that is present in the iModulonCategory column
  
  my_df <- df %>% 
    filter(iModulonCategory == Category) %>% 
    mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")"))
  
  ggplot(my_df, aes(x = LOG2FOLD, y = PathName_2)) + 
    geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                   fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                    FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                    TRUE ~ "ns")),
               size = 4, shape = 21, alpha = 0.8) + 
    scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
    guides(shape = "none") + 
    # scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) + 
    geom_vline(xintercept = 0) + 
    labs(title = paste0(deparse(substitute(df)), " ", Category), y = NULL, x = "Log2Fold change") + 
    my_plot_themes
}


###########################################################
##################### W0 Cure Vs Ra #######################

W0CureVsRa_iModulons <- MetaGeneSets_W0.cure.ComparedTo.Ra %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make new column with Groups, from Import_DEG_sets.R
W0CureVsRa_iModulons <- W0CureVsRa_iModulons %>%
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Try making one bubble plot
test <- BubblePlot_Function(W0CureVsRa_iModulons, "Central Carbon")
test

# Loop to save all categories
my_path <- paste0("Figures/Bubbles/iModulons/", "W0CureVsRa_iModulons_all")

# for (cat in unique(W0CureVsRa_iModulons$iModulonCategory)) {
#   p <- BubblePlot_Function(W0CureVsRa_iModulons, cat)
#   ggsave(p,
#          file = paste0(deparse(substitute(W0CureVsRa_iModulons)), "_", cat, "Run1-3.pdf"),
#          path = my_path,
#          width = 7.5, height = 6, units = "in")}


###########################################################
################# W0 Relapse Vs W0 Cure ###################

W0RelapseVsW0Cure_iModulons <- MetaGeneSets_W0.relapse.ComparedTo.W0.cure %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>% 
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Loop to save all categories
my_path <- paste0("Figures/Bubbles/iModulons/", "W0RelapseVsW0Cure_iModulons", "_all")

# for (cat in unique(W0RelapseVsW0Cure_iModulons$iModulonCategory)) {
#   p <- BubblePlot_Function(W0RelapseVsW0Cure_iModulons, cat)
#   ggsave(p,
#          file = paste0(deparse(substitute(W0RelapseVsW0Cure_iModulons)), "_", cat, "Run1-3.pdf"),
#          path = my_path,
#          width = 7.5, height = 6, units = "in")}


###########################################################
################### W2 Cure Vs W0 Cure ####################

W2CureVsW0Cure_iModulons <- MetaGeneSets_W2.cure.ComparedTo.W0.cure %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>% 
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Loop to save all categories
my_path <- paste0("Figures/Bubbles/iModulons/", "W2CureVsW0Cure_iModulons", "_all")

# for (cat in unique(W2CureVsW0Cure_iModulons$iModulonCategory)) {
#   p <- BubblePlot_Function(W2CureVsW0Cure_iModulons, cat)
#   ggsave(p,
#          file = paste0(deparse(substitute(W2CureVsW0Cure_iModulons)), "_", cat, "Run1-3.pdf"),
#          path = my_path,
#          width = 7.5, height = 6, units = "in")}


###########################################################
##################### W0 Relapse Vs Ra ####################

W0RelapseVsRa_iModulons <- MetaGeneSets_W0.relapse.ComparedTo.Ra %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>% 
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Loop to save all categories
my_path <- paste0("Figures/Bubbles/iModulons/", "W0RelapseVsRa_iModulons", "_all")

# for (cat in unique(W0RelapseVsRa_iModulons$iModulonCategory)) {
#   p <- BubblePlot_Function(W0RelapseVsRa_iModulons, cat)
#   ggsave(p,
#          file = paste0(deparse(substitute(W0RelapseVsRa_iModulons)), "_", cat, "Run1-3.pdf"),
#          path = my_path,
#          width = 7.5, height = 6, units = "in")}

###########################################################
###################### W2 Cure Vs Ra ######################

W2CureVsRa_iModulons <- MetaGeneSets_W2.cure.ComparedTo.Ra %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("iModulons: ", "") %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  %>%
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>% 
  mutate(iModulonCategory = case_when(
    str_detect(PathName, CentralCarbon_iModulons_pattern) ~ "Central Carbon",
    str_detect(PathName, AminoAcid_iModulons_pattern) ~ "Amino Acid",
    str_detect(PathName, NucleicAcid_iModulons_pattern) ~ "Nucleic Acid",
    str_detect(PathName, FattyAcid.Cholesterol_iModulons_pattern) ~ "Fatty Acid_Cholesterol",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, SulfurMetabolism_iModulons_pattern) ~ "Sulfur",
    str_detect(PathName, Growth_iModulons_pattern) ~ "Growth",
    str_detect(PathName, Redox_iModulons_pattern) ~ "Redox",
    str_detect(PathName, AcidStress_iModulons_pattern) ~ "Acid Stress",
    str_detect(PathName, Antibiotic_iModulons_pattern) ~ "Antibiotic",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence_Persistence",
    TRUE ~ "Other"
  ))

# Loop to save all categories
my_path <- paste0("Figures/Bubbles/iModulons/", "W2CureVsRa_iModulons", "_all")

# for (cat in unique(W2CureVsRa_iModulons$iModulonCategory)) {
#   p <- BubblePlot_Function(W2CureVsRa_iModulons, cat)
#   ggsave(p,
#          file = paste0(deparse(substitute(W2CureVsRa_iModulons)), "_", cat, "Run1-3.pdf"),
#          path = my_path,
#          width = 7.5, height = 6, units = "in")}




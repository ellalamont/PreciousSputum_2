# iModulons gene sets comparisons
# E. Lamont
# 10/13/25


# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result
# I know the UP and DOWN files are a little different, just working with the UP files for now, should maybe ask Bob if that is okay at some point

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=7),
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
facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))


###########################################################
############## CHOOSE iMODULONS OF INTEREST ###############

iModulons_of_interest <- c("PrpR", "BkaR", "Rv0681", "KstR2", "Fatty Acid Biosynthesis", "FasR", "Peptidoglycan Biosynthesis", "Positive Regulation of Growth", "Mycofactocin Synthesis Pathway", "DevR-2", "GroEL-GroES Complex", "DevR-1", "PhoP", "MprA", "Mce3R", "Mce1R", "SigC", "SigD", "SigH", "SigK", "RicR", "IdeR", "Zur", "WhiB1") 
iModulons_of_interest_pattern <- str_c(iModulons_of_interest, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or



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

# Make code simple by using my_df
my_df <- W0CureVsRa_iModulons
my_title <- "W0CureVsRa_iModulons"

Fav_Pathways <- my_df %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>%
  pull(PathName)

# Put them in my own pathway labels
my_df2 <- my_df %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR",
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))

# Make the bubble plot
my_bubblePlot <- my_df2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey")) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = my_title, y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0(my_title, ".pdf"),
#        path = "Figures/Bubbles/iModulons/fromPaper",
#        width = 6, height = 7, units = "in")


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
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make code simple by using my_df
my_df <- W0RelapseVsW0Cure_iModulons
my_title <- "W0RelapseVsW0Cure_iModulons"

Fav_Pathways <- my_df %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>%
  pull(PathName)

# Put them in my own pathway labels
my_df2 <- my_df %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR",
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))

# Make the bubble plot
my_bubblePlot <- my_df2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-2, 4), breaks = seq(-2, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = my_title, y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0(my_title, ".pdf"),
#        path = "Figures/Bubbles/iModulons/fromPaper",
#        width = 6, height = 7, units = "in")


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
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make code simple by using my_df
my_df <- W2CureVsW0Cure_iModulons
my_title <- "W2CureVsW0Cure_iModulons"

Fav_Pathways <- my_df %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>%
  pull(PathName)

# Put them in my own pathway labels
my_df2 <- my_df %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR",
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))

# Make the bubble plot
my_bubblePlot <- my_df2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = my_title, y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0(my_title, ".pdf"),
#        path = "Figures/Bubbles/iModulons/fromPaper",
#        width = 6, height = 7, units = "in")


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
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make code simple by using my_df
my_df <- W0RelapseVsRa_iModulons
my_title <- "W0RelapseVsRa_iModulons"

Fav_Pathways <- my_df %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>%
  pull(PathName)

# Put them in my own pathway labels
my_df2 <- my_df %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR",
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))

# Make the bubble plot
my_bubblePlot <- my_df2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3, 4), breaks = seq(-3, 4, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = my_title, y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0(my_title, ".pdf"),
#        path = "Figures/Bubbles/iModulons/fromPaper",
#        width = 6, height = 7, units = "in")

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
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))

# Make code simple by using my_df
my_df <- W2CureVsRa_iModulons
my_title <- "W2CureVsRa_iModulons"

Fav_Pathways <- my_df %>% 
  filter(str_detect(PathName, iModulons_of_interest_pattern)) %>%
  pull(PathName)

# Put them in my own pathway labels
my_df2 <- my_df %>%
  mutate(iModulonCategory2 = case_when(
    str_detect(PathName, "DevR") ~ "DosR",
    str_detect(PathName, paste(Growth_iModulons_pattern, NucleicAcid_iModulons_pattern, Redox_iModulons_pattern, AminoAcid_iModulons_pattern, sep = "|")) ~ "Growth",
    str_detect(PathName, Metal_iModulons_pattern) ~ "Metal",
    str_detect(PathName, Virulence.Persistence_iModulons_pattern) ~ "Virulence and Persistence",
    str_detect(PathName, paste(CentralCarbon_iModulons_pattern, FattyAcid.Cholesterol_iModulons_pattern, sep = "|")) ~ "Fatty Acid and Cholesterol",
    TRUE ~ "Other"
  )) %>%
  mutate(iModulonCategory2 = factor(iModulonCategory2, levels = c("Fatty Acid and Cholesterol", "Growth", "DosR", "Metal","Virulence and Persistence", "Other")))

# Make the bubble plot
my_bubblePlot <- my_df2 %>%
  filter(PathName %in% Fav_Pathways) %>%
  filter(iModulonCategory2 != "Virulence and Persistence") %>% 
  filter(!str_detect(PathName, "WhiB4/IdeR")) %>% # To remove this iModulon which is tagging along with IdeR iModulon
  filter(N_Genes >=3) %>% 
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.8) + 
  scale_fill_manual(
    values = c("pos" = "#bb0c00", 
               "neg" = "#00AFBB", 
               "ns"  = "grey"),
    name = "Significance / Direction"
  ) +
  facet_grid(rows = vars(iModulonCategory2), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-4, 3), breaks = seq(-4, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = my_title, y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0(my_title, ".pdf"),
#        path = "Figures/Bubbles/iModulons/fromPaper",
#        width = 6, height = 7, units = "in")




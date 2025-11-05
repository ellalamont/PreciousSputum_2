# Specific gene sets comparisons I have chosen
# E. Lamont
# 10/24/25 Run1-3


# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result
# I know the UP and DOWN files are a little different, just working with the UP files for now, should maybe ask Bob if that is okay at some point

source("Import_DEG_sets.R")
source("Import_GeneSets.R")

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10)# , 
        # plot.margin = margin(10, 10, 10, 20)# ,
  )
facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))

###########################################################
############## IMPORT GENESET GROUP METADATA ##############

EllaGeneSets_2025.10.24 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2025.10.24.csv")

###########################################################
################### W2 Cure Vs W0 Cure ####################

W2CureVsW0Cure_EllaGeneSets <- EllaGeneSets_2025.10.24_W2.cure.ComparedTo.W0.cure %>% 
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim()) %>%
  mutate(PathName = str_wrap(PathName, width = 50)) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>%
  left_join(EllaGeneSets_2025.10.24 %>% # Add the Group names
              rename(PathName = GeneSet) %>%
              select(PathName, Group), by = "PathName")

# Add the Virulence/persistence iModulons
# source("Bubble_iModulons_all.R")
# VP_iModulons_W2CureVsW0Cure <- W2CureVsW0Cure_iModulons %>%
#   filter(iModulonCategory == "Virulence_Persistence") %>%
#   rename(Group = iModulonCategory)
# W2CureVsW0Cure_EllaGeneSets <- merge(W2CureVsW0Cure_EllaGeneSets, VP_iModulons_W2CureVsW0Cure, all = T)

# Make the bubble plot
my_bubblePlot <- W2CureVsW0Cure_EllaGeneSets %>%
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>%
  mutate(Group_wrapped = case_when(Group_wrapped == "Ribosomal proteins" ~ "Ribosomal\nproteins", Group_wrapped == "Hypoxia related" ~ "Hypoxia\nrelated", TRUE ~ Group_wrapped)) %>%
  filter(Group != "Toxin/Antitoxin") %>% # Remove this because I don't think its interesting
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.9) + 
  scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
  facet_grid(rows = vars(Group_wrapped), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "W2CureVsW0Cure (Run1-3)", y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0("W2CureVsW0Cure_EllaGeneSets_2024.10.24", ".pdf"),
#        path = "Figures/Bubbles/EllaGeneSets",
#        width = 8, height = 8, units = "in")



###########################################################
################ 11/5/25 W2 Relapse vs Cure ###############

EllaGeneSets_2025.11.05 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2025.11.05.csv")











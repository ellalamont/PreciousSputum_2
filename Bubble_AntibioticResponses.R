# Compare gene sets involved in antibiotic responses
# E. Lamont
# 10/13/25
# 11/10/25 changed to be with the >60% txn cov genes run 1-3

# Gene sets were from:
# Poonawala et al. (2024)  Transcriptomic responses to antibiotic exposure in Mycobacterium tuberculosis
# https://journals.asm.org/doi/full/10.1128/aac.01185-23

# Before this was run, had to make a new gene set, add it to the lenovo, then run the MetaGeneSets function on it. Then import that back here.

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
                      strip.text = element_text(size = 10))

###########################################################
################### W2 Cure Vs W0 Cure ####################

# TAR_Poonawala2024_W2.cure.ComparedTo.W0.cure
# (TAR = transcriptional antibiotic response)

W2CureVsW0Cure_TAR <- TAR_Poonawala2024_W2.cure.ComparedTo.W0.cure %>% 
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant")) %>%
  mutate(Category = case_when(
    str_detect(PathName, "Pyrazinamide") ~ "PZA",
    str_detect(PathName, "Ethambutol") ~ "EMB",
    str_detect(PathName, "Rifampicin") ~ "RIF",
    str_detect(PathName, "Isoniazid") ~ "INH")) %>%
  mutate(Category = factor(Category, levels = c("RIF", "INH", "EMB", "PZA")))

# Make the bubble plot
my_bubblePlot <- W2CureVsW0Cure_TAR %>%
  filter(!str_detect(PathName, " 24hr")) %>% # Don't look at the >24hr set
  mutate(PathName_2 = paste0(PathName, " (n=", N_Genes, ")")) %>%
  ggplot(aes(x = LOG2FOLD, y = PathName_2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & LOG2FOLD>0 ~ "pos",
                                  FDR_Significance == "significant" & LOG2FOLD<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.9) + 
  scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
  facet_grid(rows = vars(Category), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "W2CureVsW0Cure_TAR (Run1-3)", y = NULL, x = "Log2Fold change") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
ggsave(my_bubblePlot,
       file = paste0("W2CureVsW0Cure_TAR", ".pdf"),
       path = "Figures/Bubbles/AntibioticResponse",
       width = 6, height = 5, units = "in")




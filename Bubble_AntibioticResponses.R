# Compare gene sets invovled in antibiotic responses
# E. Lamont
# 10/13/25

# Gene sets were from:
# Poonawala et al. (2024)  Transcriptomic responses to antibiotic exposure in Mycobacterium tuberculosis
# https://journals.asm.org/doi/full/10.1128/aac.01185-23

# Before this was run, had to make a new gene set, add it to the lenovo, then run the MetaGeneSets function on it. Then import that back here.

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

###########################################################
################### W2 Cure Vs W0 Cure ####################

# TAR_Poonawala2024_W2.cure.ComparedTo.W0.cure
# (TAR = transcriptional antibiotic response)

W2CureVsW0Cure_TAR <- TAR_Poonawala2024_W2.cure.ComparedTo.W0.cure %>% 
  select(PathName, CellType, N_Genes, LOG2FOLD, AVG_PVALUE, AVG_RANK) %>%
  mutate(FDR.pvalue  = p.adjust(AVG_PVALUE, method = "fdr")) %>% 
  mutate(FDR_Significance = ifelse(FDR.pvalue < 0.05, "significant", "not significant"))


# Gave up on this for now... Will try the forest plots first...






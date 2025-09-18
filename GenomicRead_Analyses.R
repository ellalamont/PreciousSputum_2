# Look at the P_Genomic and N_Genomic for the samples
# E. Lamont
# 9/18/25


source("Import_data.R") # To get All_pipeSummary


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        # axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default

###########################################################
###################### RUN 1 ONLY #########################

N_Genomic_Fig1 <- Run2_pipeSummary %>% 
  filter(Type != "NA") %>%
  ggplot(aes(x = Type, y = N_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  # geom_point(shape = 16, alpha = 0.8, size = 1.5, position = position_jitter(0.2)) + 
  # geom_hline(yintercept = 4499*0.8, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 0.9, y = 4499*0.8, label = "80%", 
           # hjust = 1.1, vjust = -0.5, color = "black") + 
  # geom_hline(yintercept = 4499*0.5, linetype = "dashed", alpha = 0.5) + 
  # annotate("text", x = 0.9, y = 4499*0.5, label = "50%", 
           # hjust = 1.1, vjust = -0.5, color = "black") + 
  geom_hline(yintercept = 1000000, linetype = "dashed", alpha = 0.5) + 
  # scale_y_continuous(limits = c(0,12000000), breaks = seq(0, 12000000, 2000000)) +
  labs(title = "PredictTB Run 2: N_Genomic for all sample types",
       x = "Sample type", 
       y = "Number of reads mapping to H37Rv") + 
  # scale_y_continuous(limits = c(0,4500), breaks = seq(0, 4500, 500)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
N_Genomic_Fig1
# ggsave(N_Genomic_Fig1,
#        file = paste0("Run2_N_Genomic_Fig1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 8, height = 5, units = "in")



###########################################################
################ P_Genomic vs SAMPLE TYPE #################

P_Genomic_Fig1 <- Run2_pipeSummary %>% 
  filter(Type != "NA") %>%
  ggplot(aes(x = Type, y = P_Genomic)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  labs(title = "PredictTB Run 2: P_Genomic for all sample types",
       x = "Sample type", 
       y = "Percent of reads mapping to H37Rv") + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
P_Genomic_Fig1
# ggsave(P_Genomic_Fig1,
#        file = paste0("Run2_P_Genomic_Fig1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 8, height = 5, units = "in")


TenReads_Fig1 <- Run2_pipeSummary %>% 
  filter(Type != "NA") %>%
  mutate(Txn_Coverage = (AtLeast.10.Reads/4499)*100) %>%
  ggplot(aes(x = Type, y = Txn_Coverage)) + 
  geom_boxplot(fill="grey", width = 0.6, outlier.size = 0.9, alpha = 0.2) + 
  geom_point(aes(fill = Type), shape = 21, alpha = 0.8, size = 1.7, position = position_jitter(0.2)) + 
  geom_hline(yintercept = 80, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) + 
  labs(title = "PredictTB Run 2: Genes with >= 10 reads aligning for all sample types",
       x = "Sample type", 
       y = "# of genes with at least 10 reads aligning") + 
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, 10)) + 
  # scale_x_continuous(trans = "log10") + 
  my_plot_themes
TenReads_Fig1
# ggsave(TenReads_Fig1,
#        file = paste0("Run2_TenReads_Fig1.pdf"),
#        path = "Figures/GenomicRead_Analyses",
#        width = 8, height = 5, units = "in")









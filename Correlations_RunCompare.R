# Compare the same sample between different runs

source("Import_data.R") 


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=14, vjust=1, hjust=1),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )


###########################################################
############ TPM_F: COLLECT DATA OF INTEREST ##############

# Genes have been filtered to keep only protein coding Rv genes and then TPM done manually (not Bob's pipeline)

# Log10 transform the data
All_tpmf_Log10 <- All_tpm_f %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Make Gene a column
All_tpmf_Log10 <- All_tpmf_Log10 %>% rownames_to_column("Gene")


###########################################################
######################## RAW READS ########################

# Genes have been filtered to keep only protein coding Rv genes

# Log10 transform the data
All_RawReadsf_Log10 <- All_RawReads_f %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Make Gene a column
All_RawReadsf_Log10 <- All_RawReadsf_Log10 %>% rename(Gene = X)


###########################################################
############### TPM_F THP1 RUN1 vs RUN2 ###################

# Using all the genes
Sample1 <- "Run1_THP1_1e6_1" # Run1
Sample2 <- "Run2_THP1_1e6_1" # Run2
ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_Compare.Run1_Run2_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############### TPM_F THP1 RUN1 vs RUN2.5 #################

# Using all the genes
Sample1 <- "Run1_THP1_1e6_1" # Run1
Sample2 <- "Run2.5_THP1_1e6_Probe20241210" # Run2
ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_Compare.Run1_Run2.5_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############## TPM_F THP1 RUN2.5 vs RUN2.5 ################

# Using all the genes
Sample1 <- "Run2.5_THP1_1e6_Probe20250912" # Run1
Sample2 <- "Run2.5_THP1_1e6_Probe20241210" # Run2
ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_Compare.Run2.5_Run2.5_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############## TPM_F THP1 RUN2.5 vs RUN3 ################

# Using all the genes
Sample1 <- "Run2.5_THP1_1e6_Probe20241210" # Run2.5
Sample2 <- "Run3_THP1_1e6" # Run3
ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_Compare.Run2.5_Run3_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
##################### TPM_F W0_12043 ######################

Sample1 <- "Run1_W0_12043" # Run1
Sample2 <- "Run2_W0_12043" # Run2
ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_12043_Compare.Run1_Run2_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
##################### TPM_F W0_14136 ######################

Sample1 <- "Run1_W0_14136" # Run1
Sample2 <- "Run2_W0_14136" # Run2

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_14136_Compare.Run1_Run2_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")




###########################################################
############### RAWREADS_F THP1 1e6 SPIKED ################

Sample1 <- "Run1_THP1_1e6_1" # Run1
Sample2 <- "Run2_THP1_1e6_1" # Run2
ScatterCorr <- All_RawReadsf_Log10 %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "RawReads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_Compare.Run1_Run2_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
################## RAWREADS_F W0_12043 ####################

Sample1 <- "Run1_W0_12043" # Run1
Sample2 <- "Run2_W0_12043" # Run2
ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_12043_Compare.Run1_Run2_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
################## RAWREADS_F W0_14136 ####################

Sample1 <- "Run1_W0_14136" # Run1
Sample2 <- "Run2_W0_14136" # Run2
ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_14136_Compare.Run1_Run2_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############### W0_13026 TPM_F + RAW READS ################

Sample1 <- "Run1_W0_13026" 
Sample2 <- "Run3_W0_13026" 

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_13026_CompareRuns_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W0_13026_CompareRuns_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")


###########################################################
############### W2_11058 TPM_F + RAW READS ################

Sample1 <- "Run3_W2_11058" 
Sample2 <- "Run2_W2_11058" 

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_11058_CompareRuns_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_11058_CompareRuns_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############### W2_12010 TPM_F + RAW READS ################

Sample1 <- "Run2_W2_12010" 
Sample2 <- "Run3_W2_12010" 

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_12010_CompareRuns_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_12010_CompareRuns_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############### W2_12029 TPM_F + RAW READS ################

Sample1 <- "Run3_W2_12029" 
Sample2 <- "Run2_W2_12029" 

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_12029_CompareRuns_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_12029_CompareRuns_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

###########################################################
############### W2_15017 TPM_F + RAW READS ################

Sample1 <- "Run3_W2_15017" 
Sample2 <- "Run2_W2_15017" 

ScatterCorr <- All_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "tpm_f; Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_15017_CompareRuns_tpmf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")

ScatterCorr <- All_RawReadsf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Raw Reads_f, Pearson correlation",
       x = paste0("Log10(RawReads+1) ", Sample1), y = paste0("Log10(RawReads+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
# ggsave(ScatterCorr,
#        file = paste0("W2_15017_CompareRuns_RawReadsf.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")



# Compare the same sample between different runs

source("Import_data.R") # To get Run1_tpm and Run2_tpm


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
################ COLLECT DATA OF INTEREST #################

# Collecting all the duplicates between the two runs
Run1_Subset <- Run1_tpm %>% select(X, THP1_1e6_1_S67, W0_12043_S32, W0_12082_S45, W0_13094_S46, W0_14136_S52, W0_15072_S48, W0_15083_S50) %>%
  rename(Gene = X, Run1_THP1.Ra1e6 = THP1_1e6_1_S67, Run1_W0_12043 = W0_12043_S32, Run1_W0_12082 = W0_12082_S45,Run1_W0_13094 = W0_13094_S46,  Run1_W0_14136 = W0_14136_S52, Run1_W0_15072 = W0_15072_S48, Run1_W0_15083 = W0_15083_S50)
Run2_Subset <- Run2_tpm %>% select(X, THP1_1e6_1_S44, W0_12043_S45, W0_12082_S46, W0_13094_S47, W0_14136_S50, W0_15072_S48, W0_15083_S49) %>%
  rename(Gene = X, Run2_THP1.Ra1e6 = THP1_1e6_1_S44, Run2_W0_12043 = W0_12043_S45, Run2_W0_12082 = W0_12082_S46, Run2_W0_13094 = W0_13094_S47, Run2_W0_14136 = W0_14136_S50, Run2_W0_15072 = W0_15072_S48, Run2_W0_15083 = W0_15083_S49)

merged_DoubleRun <- merge(Run1_Subset, Run2_Subset, all = T)

# Log10 transform the data
merged_DoubleRun_Log10 <- merged_DoubleRun %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Remove all the non Rv genes and see how much better it gets
merged_DoubleRun_Log10_filtered <- merged_DoubleRun_Log10 %>% 
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", Gene))

###########################################################
############# GGCORRPLOT OF THP1 1e6 SPIKED ###############

# Using all the genes
Sample1 <- "Run1_THP1.Ra1e6" # Run1
Sample2 <- "Run2_THP1.Ra1e6" # Run2
ScatterCorr <- merged_DoubleRun_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("THP1 spiked with 1e6 on two separate runs (different mRNA+LibraryPrep)"),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr
ggplotly(ScatterCorr)
# ggsave(ScatterCorr,
#        file = paste0("THP1Spiked1e6_CompareRuns_AllGenes.pdf"),
#        path = "Figures/Correlations_RunCompare",
#        width = 7, height = 5, units = "in")



###########################################################
######################## W0_12043 #########################

Sample1 <- "Run1_W0_12043" # Run1
Sample2 <- "Run2_W0_12043" # Run2
ScatterCorr <- merged_DoubleRun_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr

# With just the Rv genes
Sample1 <- "Run1_W0_12043" # Run1
Sample2 <- "Run2_W0_12043" # Run2
ScatterCorr_f <- merged_DoubleRun_Log10_filtered %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Only coding Rv genes shown, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr_f
ggplotly(ScatterCorr_f)


###########################################################
######################## W0_14136 #########################

Sample1 <- "Run1_W0_14136" # Run1
Sample2 <- "Run2_W0_14136" # Run2

ScatterCorr <- merged_DoubleRun_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr

# With just the Rv genes
ScatterCorr_f <- merged_DoubleRun_Log10_filtered %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("Sputum sample sequenced on two separate runs (different capture)"),
       subtitle = "Only coding Rv genes shown, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr_f



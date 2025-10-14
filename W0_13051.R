# Look at the sequencing for the 2 (will be 3) times I sequenced W0_13051 and see what's going on with the high N_Genomic bu low txn coverage
# E. Lamont
# 10/14/25

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

# Stop scientific notation
options(scipen = 999) 


W0_13051_pipeSummary <- All_pipeSummary %>% filter(Patient == "P_13051")
W0_13051_tpmf <- All_tpm_f %>% select(contains("13051"))
W0_13051_RawReadsf <- All_RawReads_f %>% select(X, contains("13051"))


###########################################################
############### NUMBER OF GENES WITH 9 READS ##############

Num1Reads <- colSums(W0_13051_RawReadsf >= 1)
# X   Run2_W0_13051 Run2.5_W0_13051 
# 4030            3485            3702 

Num10Reads <- colSums(W0_13051_RawReadsf >= 10)
# X   Run2_W0_13051 Run2.5_W0_13051 
# 4030            1861            1897 

###########################################################
#################### TPM_F CORRELATION ####################

# Genes have been filtered to keep only protein coding Rv genes and then TPM done manually (not Bob's pipeline)

# Log10 transform the data
W0_13051_tpmf_Log10 <- W0_13051_tpmf %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Make Gene a column
W0_13051_tpmf_Log10 <- W0_13051_tpmf_Log10 %>% rownames_to_column("Gene")

Sample1 <- "Run2_W0_13051" # Run1
Sample2 <- "Run2.5_W0_13051" # Run2
ScatterCorr <- W0_13051_tpmf_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  # geom_text_repel(aes(label = Gene), size= 0.5, max.overlaps = 3) + 
  geom_text(aes(label = Gene), size = 2, vjust = -0.5, hjust = 0.5, check_overlap = T) +  
  labs(title = paste0("W0_13051 sequenced twice"),
       subtitle = "tpm_f, Pearson correlation",
       x = paste0("Log10(TPM+1) ", Sample1), y = paste0("Log10(TPM+1) ", Sample2)) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  # scale_x_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  # scale_y_continuous(limits = c(0,14000), breaks = seq(0, 14000, 2000)) + 
  my_plot_themes
ScatterCorr




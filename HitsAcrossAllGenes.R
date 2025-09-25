# See what genes the reads are falling on and if they are the same across the same sample run twice
# E. Lamont
# 9/23/25

source("Import_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
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
############## TPM: COLLECT DATA OF INTEREST ##############

# Collecting all the duplicates between the two runs
tpm_Run1_Subset <- Run1_tpm %>% select(X, THP1_1e6_1_S67, W0_12043_S32, W0_12082_S45, W0_13094_S46, W0_14136_S52, W0_15072_S48, W0_15083_S50) %>%
  rename(Gene = X, Run1_THP1.Ra1e6 = THP1_1e6_1_S67, Run1_W0_12043 = W0_12043_S32, Run1_W0_12082 = W0_12082_S45,Run1_W0_13094 = W0_13094_S46,  Run1_W0_14136 = W0_14136_S52, Run1_W0_15072 = W0_15072_S48, Run1_W0_15083 = W0_15083_S50)
tpm_Run2_Subset <- Run2_tpm %>% select(X, THP1_1e6_1_S44, W0_12043_S45, W0_12082_S46, W0_13094_S47, W0_14136_S50, W0_15072_S48, W0_15083_S49) %>%
  rename(Gene = X, Run2_THP1.Ra1e6 = THP1_1e6_1_S44, Run2_W0_12043 = W0_12043_S45, Run2_W0_12082 = W0_12082_S46, Run2_W0_13094 = W0_13094_S47, Run2_W0_14136 = W0_14136_S50, Run2_W0_15072 = W0_15072_S48, Run2_W0_15083 = W0_15083_S49)

tpm_merged_DoubleRun <- merge(tpm_Run1_Subset, tpm_Run2_Subset, all = T)

# Log10 transform the data
tpm_DoubleRun_Log10 <- tpm_merged_DoubleRun %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Pivot longer the data
tpm_DoubleRun_Log10_t <- tpm_DoubleRun_Log10 %>%
  pivot_longer(cols = starts_with("Run"), names_to = "Sample", values_to = "TPM_log10")



###########################################################
########### RAW READS: COLLECT DATA OF INTEREST ###########

# Collecting all the duplicates between the two runs
RawReads_Run1_Subset <- Run1_RawReads %>% select(X, THP1_1e6_1_S67, W0_12043_S32, W0_12082_S45, W0_13094_S46, W0_14136_S52, W0_15072_S48, W0_15083_S50) %>%
  rename(Gene = X, Run1_THP1.Ra1e6 = THP1_1e6_1_S67, Run1_W0_12043 = W0_12043_S32, Run1_W0_12082 = W0_12082_S45,Run1_W0_13094 = W0_13094_S46,  Run1_W0_14136 = W0_14136_S52, Run1_W0_15072 = W0_15072_S48, Run1_W0_15083 = W0_15083_S50)
RawReads_Run2_Subset <- Run2_RawReads %>% select(X, THP1_1e6_1_S44, W0_12043_S45, W0_12082_S46, W0_13094_S47, W0_14136_S50, W0_15072_S48, W0_15083_S49) %>%
  rename(Gene = X, Run2_THP1.Ra1e6 = THP1_1e6_1_S44, Run2_W0_12043 = W0_12043_S45, Run2_W0_12082 = W0_12082_S46, Run2_W0_13094 = W0_13094_S47, Run2_W0_14136 = W0_14136_S50, Run2_W0_15072 = W0_15072_S48, Run2_W0_15083 = W0_15083_S49)

RawReads_merged_DoubleRun <- merge(RawReads_Run1_Subset, RawReads_Run2_Subset, all = T)

# Remove all the non Rv genes and see how much better it gets
RawReads_DoubleRun_filtered <- RawReads_merged_DoubleRun %>% 
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", Gene))

# Pivot longer the data
RawReads_DoubleRun_ft <- RawReads_DoubleRun_filtered %>%
  pivot_longer(cols = starts_with("Run"), names_to = "Sample", values_to = "RawReads")

###########################################################
######################### THP1 SPIKED #####################

# TPM
THP1_1e6_fig <- tpm_DoubleRun_Log10_t %>%
  filter(str_detect(Sample, "THP1")) %>% 
  # head(20) %>% 
  ggplot(aes(x = Gene, y = TPM_log10, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Genes hits for THP1 1e6 spiked samples",
       # subtitle = "big gene is Rv2503c", 
       y = "log10(TPM + 1)",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes
# THP1_1e6_fig # too big to print in R
# ggplotly(THP1_1e6_fig)
# ggsave(THP1_1e6_fig,
#        file = paste0("THP1_1e6_fig.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 20, height = 5, units = "in")

# Filtered Raw Reads
THP1_1e6_RawReads_v1 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "THP1")) %>% 
  # head(20) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Genes (Rv only) hits for THP1 1e6 spiked samples",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes 
# ggsave(THP1_1e6_RawReads_v1,
#        file = paste0("THP1_1e6_RawReads_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 20, height = 5, units = "in")

# Top 100 genes
top100Genes_Run1 <- RawReads_DoubleRun_filtered %>% 
  select(Gene, Run1_THP1.Ra1e6) %>%
  arrange(desc(Run1_THP1.Ra1e6)) %>%
  slice_head(n=100) %>%
  pull(Gene)

THP1_1e6_RawReads_top100_v1 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "THP1")) %>% 
  filter(Gene %in% top100Genes_Run1) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Top 100 most expressed Rv genes for THP1 1e6 spiked samples",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes + theme(axis.text.x = element_text(angle = 90, size=8, vjust=0.5, hjust=1))
THP1_1e6_RawReads_top100_v1
# ggsave(THP1_1e6_RawReads_top100_v1,
#        file = paste0("THP1_1e6_RawReads_top100_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 15, height = 5, units = "in")


###########################################################
######################### W0_12043 ########################

W0_12043_fig <- tpm_DoubleRun_Log10_t %>%
  filter(str_detect(Sample, "W0_12043")) %>% 
  # head(20) %>% 
  ggplot(aes(x = Gene, y = TPM_log10, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Genes hits for W0_12043",
       # subtitle = "big gene is Rv2503c", 
       y = "log10(TPM + 1)",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes
# ggsave(W0_12043_fig,
#        file = paste0("W0_12043_fig.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 20, height = 5, units = "in")


# Filtered Raw Reads
W0_12043_RawReads_v1 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "W0_12043")) %>% 
  # head(20) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Genes (Rv only) hits for W0_12043",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes 
# ggsave(W0_12043_RawReads_v1,
#        file = paste0("W0_12043_RawReads_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 20, height = 5, units = "in")

# Top 100 genes
top100Genes_Run1 <- RawReads_DoubleRun_filtered %>% 
  select(Gene, Run1_W0_12043) %>%
  arrange(desc(Run1_W0_12043)) %>%
  slice_head(n=100) %>%
  pull(Gene)

W0_12043_RawReads_top100_v1 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "W0_12043")) %>% 
  filter(Gene %in% top100Genes_Run1) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Top 100 most expressed Rv genes for W0_12043",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes + theme(axis.text.x = element_text(angle = 90, size=8, vjust=0.5, hjust=1))
W0_12043_RawReads_top100_v1
# ggsave(THP1_1e6_RawReads_top100_v1,
#        file = paste0("THP1_1e6_RawReads_top100_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 15, height = 5, units = "in")


###########################################################
######################### W0_12082 ########################


# Top 50 genes Run1
top50Genes_Run1 <- RawReads_DoubleRun_filtered %>% 
  select(Gene, Run1_W0_12082) %>%
  arrange(desc(Run1_W0_12082)) %>%
  slice_head(n=50) %>%
  pull(Gene)

W0_12082_RawReads_top50_byRun1 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "W0_12082")) %>% 
  filter(Gene %in% top50Genes_Run1) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Top 50 most expressed Rv genes in Run1 for W0_12082",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes + theme(axis.text.x = element_text(angle = 90, size=8, vjust=0.5, hjust=1))
W0_12082_RawReads_top50_byRun1
# ggsave(THP1_1e6_RawReads_top100_v1,
#        file = paste0("THP1_1e6_RawReads_top100_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 15, height = 5, units = "in")

# Top 50 genes Run2
top50Genes_Run2 <- RawReads_DoubleRun_filtered %>% 
  select(Gene, Run2_W0_12082) %>%
  arrange(desc(Run2_W0_12082)) %>%
  slice_head(n=50) %>%
  pull(Gene)

W0_12082_RawReads_top50_byRun2 <- RawReads_DoubleRun_ft %>%
  filter(str_detect(Sample, "W0_12082")) %>% 
  filter(Gene %in% top50Genes_Run2) %>% 
  ggplot(aes(x = Gene, y = RawReads, fill = Sample)) + # for log10 graph
  # ggplot(aes(x = X, y = TPM, fill = Sample), text = X) + # for normal graph
  geom_bar(position="fill", stat="identity") +  # position = stack for raw reads, position = fill for propotion
  scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
  labs(title = "Top 50 most expressed Rv genes in Run2 for W0_12082",
       # subtitle = "big gene is Rv2503c", 
       y = "RawReads",
       x = "Gene (alphabetically ordered)") + 
  scale_y_continuous(expand = c(0,0)) + 
  my_plot_themes + theme(axis.text.x = element_text(angle = 90, size=8, vjust=0.5, hjust=1))
W0_12082_RawReads_top50_byRun2
# ggsave(THP1_1e6_RawReads_top100_v1,
#        file = paste0("THP1_1e6_RawReads_top100_v1.pdf"),
#        path = "Figures/HitsAcrossAllGenes",
#        width = 15, height = 5, units = "in")

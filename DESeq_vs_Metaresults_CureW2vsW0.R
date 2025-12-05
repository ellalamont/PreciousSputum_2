# Testing the difference between DESeq2 and metaresults for Cure W0 and W2
# 12/2/25

source("Import_data.R")
source("Import_GeneSets.R")

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering


################################################
################ LOAD PACKAGES #################

# https://github.com/JulioLeonIncio/Tutorial-from-DEseq2-to-GSEA-in-R/blob/main/DEseq2_to_GSEA_JL.Rmd
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# pkgs <- c("DESeq2","sva","fgsea","clusterProfiler","GSEABase","tidyverse","pheatmap", "apeglm")
# for (p in pkgs) {
#   if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE)
# }
library(DESeq2)
library(sva) # ComBat_seq
library(fgsea)
library(clusterProfiler)
library(GSEABase)
library(tidyverse)
library(pheatmap)
library(fgsea)
library(dplyr)
library(msigdbr)
library(apeglm)


################################################
################## PRE-FILTER ##################
# Remove genes with <10 TPM

# sputum_tpm <- GoodSamples60_tpmf %>% dplyr::select(contains("W"))
# goodGenesList <- sputum_tpm
# 
# goodGenesList <- rownames(GoodSamples60_tpmf)[apply(GoodSamples60_tpmf, 1, function(x) all(x >= 10))]

################################################
##################### DESEQ2 ###################

# Keeping just the sputum samples at the beginning
# Bob's code: https://github.com/robertdouglasmorrison/DuffyTools/blob/master/R/DESeqTools.R
# Using local fit (because this seems to be what bob's code is doing)
# MetaResults has cooksCutoff set to FALSE, I think this may fix the NA p-values
##### Yes, setting it to FALSE, removed the NA's from the p-values
# Also set independentFiltering = F
# Not clear on how to shrink the log2fold values


# Start with GoodSamples60_RawReadsf
my_RawCounts <- GoodSamples60_RawReadsf %>% # rename("GENE_ID" = "X")
  column_to_rownames("X") %>%
  dplyr::select(contains("W")) %>%
  mutate_if(is.numeric, ~ as.integer(round(., 0)))

# Start with GoodSamples60_pipeSummary
my_metadata <- GoodSamples60_pipeSummary %>% 
  filter(Type %in% c("Week 0 sputum", "Week 2 sputum")) %>% # Keep just the sputum samples
  dplyr::select(SampleID, SampleID2, Type2) %>%
  rename("sample" = "SampleID2", "condition" = "Type2") %>%
  column_to_rownames("sample")
my_metadata$condition <- gsub(" ", "_", my_metadata$condition)
my_metadata$condition <- gsub("\\(|\\)", "", my_metadata$condition)
my_metadata$condition <- factor(my_metadata$condition)

# Ensure sample names line up
stopifnot(all(colnames(my_RawCounts) == rownames(my_metadata)))
stopifnot(all(rownames(my_metadata) == colnames(my_RawCounts)))

# Build DESeq
dds <- DESeqDataSetFromMatrix(countData = my_RawCounts,
                              colData = my_metadata,
                              design = ~ condition)

# Make W0_sputum_cure the reference (need to do this to shrink the log2fold changes)
dds$condition <- relevel(dds$condition, ref = "W0_sputum_cure")

dds <- DESeq(dds, fitType = "local") 

resultsNames(dds) # see what my coefficients are 

# Shrink LFC to mimic DuffyTools
# res_W2Relapse_vs_W2Cure <- lfcShrink(dds, 
#                  coef = "condition_W2_sputum_relapse_vs_W2_sputum_cure",
#                  type = "apeglm") 
# summary(res_W2Relapse_vs_W2Cure)

res_W2Cure_vs_W0Cure <- results(dds, contrast = c("condition", "W2_sputum_cure", "W0_sputum_cure"), cooksCutoff = F, independentFiltering = F)
summary(res_W2Cure_vs_W0Cure)

# Save output
res_W2Cure_vs_W0Cure_df <- as.data.frame(res_W2Cure_vs_W0Cure) %>% rownames_to_column("gene")


# Columns for Log2Fold > 2
res_W2Cure_vs_W0Cure_df$DE2 <- ifelse(res_W2Cure_vs_W0Cure_df$log2FoldChange < -2 & res_W2Cure_vs_W0Cure_df$padj < 0.05, "significant down", ifelse(res_W2Cure_vs_W0Cure_df$log2FoldChange > 2 & res_W2Cure_vs_W0Cure_df$padj < 0.05, "significant up", "not significant"))
res_W2Cure_vs_W0Cure_df$DE2 <- factor(res_W2Cure_vs_W0Cure_df$DE2, levels = c("significant down", "not significant", "significant up"))
res_W2Cure_vs_W0Cure_df$DE2_labels <- ifelse(res_W2Cure_vs_W0Cure_df$DE2 != "not significant", res_W2Cure_vs_W0Cure_df$gene, NA)

# write.csv(res_W2Relapse_vs_W2Cure_df, "Data/DESeq2_GSEA_Output/res_W2Relapse_vs_W2Cure_df4.csv")

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=10))

my_volcano <- res_W2Cure_vs_W0Cure_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE2, label = DE2_labels, text = gene)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "DESeq2 W2 cure vs W0 cure Log2Fold=2 >60%TxnCov") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- res_W2Cure_vs_W0Cure_df %>% filter(DE2 == "significant up") %>% nrow()
text_down <- res_W2Cure_vs_W0Cure_df %>% filter(DE2 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano
# ggsave(final_volcano,
#        file = "W2.cure.ComparedTo.W0.cure_DESeq2_Run1to3_v1.pdf",
#        path = "Figures/Volcano_60TxnCov/DESeq2_ByHand",
#        width = 7, height = 5, units = "in")

# To see which genes aren't being plotted
res_W2Cure_vs_W0Cure_df %>% 
  filter(is.na(log2FoldChange) | is.na(padj) | is.na(DE2))

################################################
######### METARESULTS DESeq FDR-PVALUE #########

# Taking from the DESeq ratio file

# Do a p-value correction
`DESeq_W2.cure.ComparedTo.W0.cure` <- read.delim("Data/DE_Run1to3_60TxnCov/W0.cure_vs_W2.cure/W2_cure.MTb.DESeq.Ratio.txt") %>% 
  filter(str_detect(GENE_ID, "^Rv\\d+.*"))
DESeq_W2.cure.ComparedTo.W0.cure$FDR_PVALUE <- p.adjust(DESeq_W2.cure.ComparedTo.W0.cure$PVALUE, method = "fdr")

# Columns for Log2Fold > 2
DESeq_W2.cure.ComparedTo.W0.cure$DE2 <- ifelse(DESeq_W2.cure.ComparedTo.W0.cure$LOG2FOLD < -2 & DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE < 0.05, "significant down", ifelse(DESeq_W2.cure.ComparedTo.W0.cure$LOG2FOLD > 2 & DESeq_W2.cure.ComparedTo.W0.cure$FDR_PVALUE < 0.05, "significant up", "not significant"))
DESeq_W2.cure.ComparedTo.W0.cure$DE2 <- factor(DESeq_W2.cure.ComparedTo.W0.cure$DE2, levels = c("significant down", "not significant", "significant up"))
DESeq_W2.cure.ComparedTo.W0.cure$DE2_labels <- ifelse(DESeq_W2.cure.ComparedTo.W0.cure$DE2 != "not significant", DESeq_W2.cure.ComparedTo.W0.cure$GENE_NAME, NA)

my_volcano <- DESeq_W2.cure.ComparedTo.W0.cure %>%
  ggplot(aes(x = LOG2FOLD, y = -log10(FDR_PVALUE), col = DE2, label = DE2_labels, text = GENE_NAME)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "METARESULTS DESeq2 W2 cure vs W0 cure Log2Fold=2 >60%TxnCov") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- DESeq_W2.cure.ComparedTo.W0.cure %>% filter(DE2 == "significant up") %>% nrow()
text_down <- DESeq_W2.cure.ComparedTo.W0.cure %>% filter(DE2 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano
# ggsave(final_volcano,
#        file = "W2.cure.ComparedTo.W0.Cure_BobDESeq2_Run1to3_v1.pdf",
#        path = "Figures/Volcano_60TxnCov/Log2Fold2_FDR",
#        width = 7, height = 5, units = "in")


###################################################################
##########  COMPARE LOG2FOLD CHANGES (ADJUSTED PVALUES)  ##########

df_DESeq_inR <- res_W2Cure_vs_W0Cure_df %>% 
  dplyr::select(gene, log2FoldChange, DE2, padj) %>% 
  rename(LOG2FOLD_DESeq_inR = log2FoldChange,
         DE2_DESeq_inR = DE2,
         GENE_ID = gene,
         padj_DESeq_inR = padj) 
df_DESeq_Bob <- DESeq_W2.cure.ComparedTo.W0.cure %>% 
  dplyr::select(GENE_ID, LOG2FOLD, DE2, PVALUE, FDR_PVALUE) %>% 
  rename(LOG2FOLD_DESeq_Bob = LOG2FOLD,
         DE2_DESeq_Bob = DE2,
         PVALUE_Bob = PVALUE,
         FDR_PVALUE_Bob = FDR_PVALUE) 
df_combined <- merge(df_DESeq_inR, df_DESeq_Bob)

# Want just to show significance with colors, not necessarily log2fold change
df_combined <- df_combined %>%
  mutate(Combined_Sig = case_when(
    padj_DESeq_inR>=0.05 & FDR_PVALUE_Bob>=0.05 ~ "both not significant",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR>0 & FDR_PVALUE_Bob<0.05 & LOG2FOLD_DESeq_Bob<0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR<0 & FDR_PVALUE_Bob<0.05 & LOG2FOLD_DESeq_Bob>0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & FDR_PVALUE_Bob<0.05 ~ "both significant the same",
    padj_DESeq_inR<0.05 & FDR_PVALUE_Bob>=0.05 ~ "only one significant",
    padj_DESeq_inR>=0.05 & FDR_PVALUE_Bob<0.05 ~ "only one significant"
  ))

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        legend.box.background = element_blank())

Sample1 <- "LOG2FOLD_DESeq_inR" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_DESeq_Bob" # Broth Not Captured
my_plot <- df_combined %>% drop_na() %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID, color = Combined_Sig), alpha = 0.7, size = 2) +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("both not significant" = "#999990", "both significant the same" = "blue1", "only one significant" = "maroon4", "significant in opposite directions" = "orange")) + 
  labs(title = "W2 Cure vs W0 Cure DEG Log2Fold=2 >60%TxnCov",
       subtitle = "Pearson correlation, all p-values adjusted",
       x = Sample1, 
       y = Sample2) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot
# ggsave(my_plot,
#        file = paste0("Cure_W2vsW0_AdjP_v1.pdf"),
#        path = "Figures/DEG_Comparisons_60TxnCov",
#        width = 9, height = 6, units = "in")

# Summarize the number of genes in each group
df_combined %>%
  group_by(Combined_Sig) %>%
  summarize(n = n())



###################################################################
##########  COMPARE LOG2FOLD CHANGES (ORIGINAL PVALUES)  ##########

df_DESeq_inR <- res_W2Cure_vs_W0Cure_df %>% 
  dplyr::select(gene, log2FoldChange, DE2, padj) %>% 
  rename(LOG2FOLD_DESeq_inR = log2FoldChange,
         DE2_DESeq_inR = DE2,
         GENE_ID = gene,
         padj_DESeq_inR = padj) 
df_DESeq_Bob <- DESeq_W2.cure.ComparedTo.W0.cure %>% 
  dplyr::select(GENE_ID, LOG2FOLD, DE2, PVALUE, FDR_PVALUE) %>% 
  rename(LOG2FOLD_DESeq_Bob = LOG2FOLD,
         DE2_DESeq_Bob = DE2,
         PVALUE_Bob = PVALUE,
         FDR_PVALUE_Bob = FDR_PVALUE) 
df_combined <- merge(df_DESeq_inR, df_DESeq_Bob)

# Want just to show significance with colors, not necessarily log2fold change
# Changing this to be the unadjusted p-values from Bob's pipeline
df_combined <- df_combined %>%
  mutate(Combined_Sig = case_when(
    padj_DESeq_inR>=0.05 & PVALUE_Bob>=0.05 ~ "both not significant",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR>0 & PVALUE_Bob<0.05 & LOG2FOLD_DESeq_Bob<0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR<0 & PVALUE_Bob<0.05 & LOG2FOLD_DESeq_Bob>0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & PVALUE_Bob<0.05 ~ "both significant the same",
    padj_DESeq_inR<0.05 & PVALUE_Bob>=0.05 ~ "only one significant",
    padj_DESeq_inR>=0.05 & PVALUE_Bob<0.05 ~ "only one significant"
  ))

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        legend.box.background = element_blank())

Sample1 <- "LOG2FOLD_DESeq_inR" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_DESeq_Bob" # Broth Not Captured
my_plot <- df_combined %>% drop_na() %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID, color = Combined_Sig), alpha = 0.7, size = 2) +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("both not significant" = "#999990", "both significant the same" = "blue1", "only one significant" = "maroon4", "significant in opposite directions" = "orange")) + 
  labs(title = "W2 Cure vs W0 Cure DEG Log2Fold=2 >60%TxnCov",
       subtitle = "Pearson correlation, Bob's p-value NOT adjusted",
       x = Sample1, 
       y = Sample2) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot
# ggsave(my_plot,
#        file = paste0("Cure_W2vsW0_ogP_v1.pdf"),
#        path = "Figures/DEG_Comparisons_60TxnCov",
#        width = 9, height = 6, units = "in")
# Summarize the number of genes in each group
df_combined %>%
  group_by(Combined_Sig) %>%
  summarize(n = n())



#######################################################################
######  METARESULTS COMPARE LOG2FOLD CHANGES (ORIGINAL PVALUES)  ######
# Instead of Bob's DESeq only, use the entire metaresults

df_DESeq_inR <- res_W2Cure_vs_W0Cure_df %>% 
  dplyr::select(gene, log2FoldChange, DE2, padj) %>% 
  rename(LOG2FOLD_DESeq_inR = log2FoldChange,
         DE2_DESeq_inR = DE2,
         GENE_ID = gene,
         padj_DESeq_inR = padj) 
df_metaresults <- list_dfs_f2$W2.cure.ComparedTo.W0.cure %>% 
  dplyr::select(GENE_ID, LOG2FOLD, DE2, AVG_PVALUE, FDR_PVALUE) %>% 
  rename(LOG2FOLD_metaresults = LOG2FOLD,
         DE2_metaresults = DE2,
         PVALUE_metaresults = AVG_PVALUE,
         FDR_PVALUE_metaresults = FDR_PVALUE) 
df_combined <- merge(df_DESeq_inR, df_metaresults)

# Want just to show significance with colors, not necessarily log2fold change
# Changing this to be the unadjusted p-values from Bob's pipeline
df_combined <- df_combined %>%
  mutate(Combined_Sig = case_when(
    padj_DESeq_inR>=0.05 & PVALUE_metaresults>=0.05 ~ "both not significant",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR>0 & PVALUE_metaresults<0.05 & LOG2FOLD_metaresults<0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & LOG2FOLD_DESeq_inR<0 & PVALUE_metaresults<0.05 & LOG2FOLD_metaresults>0 ~ "significant in opposite directions",
    padj_DESeq_inR<0.05 & PVALUE_metaresults<0.05 ~ "both significant the same",
    padj_DESeq_inR<0.05 & PVALUE_metaresults>=0.05 ~ "only one significant",
    padj_DESeq_inR>=0.05 & PVALUE_metaresults<0.05 ~ "only one significant"
  ))

my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        legend.box.background = element_blank())

Sample1 <- "LOG2FOLD_DESeq_inR" # THP1 spiked Captured
Sample2 <- "LOG2FOLD_metaresults" # Broth Not Captured
my_plot <- df_combined %>% drop_na() %>%
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = GENE_ID, color = Combined_Sig), alpha = 0.7, size = 2) +
  # geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("both not significant" = "#999990", "both significant the same" = "blue1", "only one significant" = "maroon4", "significant in opposite directions" = "orange")) + 
  labs(title = "W2 Cure vs W0 Cure DEG Log2Fold=2 >60%TxnCov",
       subtitle = "Pearson correlation, Bob's p-value NOT adjusted",
       x = Sample1, 
       y = Sample2) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
my_plot
# ggsave(my_plot,
#        file = paste0("Cure_W2vsW0_MetaResults_ogP_v1.pdf"),
#        path = "Figures/DEG_Comparisons_60TxnCov",
#        width = 9, height = 6, units = "in")

# Summarize the number of genes in each group
df_combined %>%
  group_by(Combined_Sig) %>%
  summarize(n = n())




##############################################################
##################### DESeq in R to GSEA #####################
# Remove NAs 
res <- res_W2Cure_vs_W0Cure_df %>%
  filter(!is.na(padj)) %>% 
  arrange(desc(log2FoldChange))

# Create a ranked vector (using stat or Log2FoldChange)
ranks <- res$log2FoldChange
names(ranks) <- res$gene
ranks <- sort(ranks, decreasing = TRUE)

# Filter out NAs or duplicated genes
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[!duplicated(names(ranks))]

# Load custom gene sets. Needs to have columns Gene, GeneSet. Will take from allGeneSetList
custom_list <- allGeneSetList[["EllaGeneSets_2025.11.05"]]
EllaGeneSets_2025.11.05 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2025.11.05.csv") # Also need the csv to get the grouping for plotting with facets

# Run fgsea
fgsea_res <- fgsea(
  pathways = custom_list,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

###########################################################
################### DESeq in R DOT PLOT ###################
my_plot_themes <- theme_bw() +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 10),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        # axis.text.x = element_text(angle = 45, size=7, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10)
  )
facet_themes <- theme(strip.background=element_rect(fill="white", linewidth = 0.9),
                      strip.text = element_text(size = 7))

# Add groupings
fgsea_res2 <- fgsea_res %>%
  mutate(pathway = str_wrap(pathway, width = 50)) %>% 
  mutate(FDR_Significance = ifelse(padj < 0.05, "significant", "not significant")) %>%
  left_join(EllaGeneSets_2025.11.05 %>% # Add the Group names
              rename(pathway = GeneSet) %>%
              dplyr::select(pathway, Group), by = "pathway")

# Make the bubble plot
my_bubblePlot <- fgsea_res2 %>%
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>%
  mutate(Group_wrapped = case_when(Group_wrapped == "Ribosomal proteins" ~ "Ribosomal\nproteins", Group_wrapped == "Hypoxia related" ~ "Hypoxia\nrelated", TRUE ~ Group_wrapped)) %>%
  filter(!Group %in% c("Toxin/Antitoxin", "ESX genes", "Metal", "Nucleic Acid")) %>% # Remove this because I don't think its interesting
  mutate(pathway2 = paste0(pathway, " (n=", size, ")")) %>%
  ggplot(aes(x = NES, y = pathway2)) + 
  geom_point(aes(stroke = ifelse(FDR_Significance == "significant", 0.8, 0),
                 fill = case_when(FDR_Significance == "significant" & NES>0 ~ "pos",
                                  FDR_Significance == "significant" & NES<0 ~ "neg",
                                  TRUE ~ "ns")),
             size = 4, shape = 21, alpha = 0.9) + 
  scale_fill_manual(values = c("pos" = "#bb0c00", "neg" = "#00AFBB", "ns"  = "grey")) +
  facet_grid(rows = vars(Group_wrapped), scales = "free_y", space = "free") + 
  guides(shape = "none") + 
  scale_x_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, 1)) + 
  geom_vline(xintercept = 0) + 
  labs(title = "W2CureVsW0Cure (Run1-3) >60%TxnCov DESeq2->GSEA", y = NULL, x = "Normalized Enrichment Score") + 
  my_plot_themes + facet_themes + theme(legend.position = "none")
my_bubblePlot
# ggsave(my_bubblePlot,
#        file = paste0("W2CureVsW0Cure_EllaGeneSets_2025.11.05_DESeq2_GSEA_v2", ".pdf"),
#        path = "Figures/Bubbles/DESeq2_GSEA/EllaGeneSets",
#        width = 8, height = 8, units = "in")





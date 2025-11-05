# Testing the difference between DESeq2 and metaresults for W2 cure vs W2 relapse
# 11/5/25

source("Import_data.R")
source("Import_GeneSets.R")

# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering


################################################
################ LOAD PACKAGES #################

# https://github.com/JulioLeonIncio/Tutorial-from-DEseq2-to-GSEA-in-R/blob/main/DEseq2_to_GSEA_JL.Rmd
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("DESeq2","sva","fgsea","clusterProfiler","GSEABase","tidyverse","pheatmap", "apeglm")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE)
}
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
#################### DESEQ2 ####################

# Keeping just the sputum samples at the beginning
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

# Make W2_sputum_cure the reference (need to do this to shrink the log2fold changes)
dds$condition <- relevel(dds$condition, ref = "W2_sputum_cure")

dds <- DESeq(dds, fitType = "local") 

resultsNames(dds) # see what my coefficients are 

# Shrink LFC to mimic DuffyTools
# res_W2Relapse_vs_W2Cure <- lfcShrink(dds, 
#                  coef = "condition_W2_sputum_relapse_vs_W2_sputum_cure",
#                  type = "apeglm") 
# summary(res_W2Relapse_vs_W2Cure)

res_W2Relapse_vs_W2Cure <- results(dds, contrast = c("condition", "W2_sputum_relapse", "W2_sputum_cure"), cooksCutoff = F, independentFiltering = F)
summary(res_W2Relapse_vs_W2Cure)

# Save output
res_W2Relapse_vs_W2Cure_df <- as.data.frame(res_W2Relapse_vs_W2Cure) %>% rownames_to_column("gene")


# Columns for Log2Fold > 2
res_W2Relapse_vs_W2Cure_df$DE2 <- ifelse(res_W2Relapse_vs_W2Cure_df$log2FoldChange < -2 & res_W2Relapse_vs_W2Cure_df$padj < 0.05, "significant down", ifelse(res_W2Relapse_vs_W2Cure_df$log2FoldChange > 2 & res_W2Relapse_vs_W2Cure_df$padj < 0.05, "significant up", "not significant"))
res_W2Relapse_vs_W2Cure_df$DE2 <- factor(res_W2Relapse_vs_W2Cure_df$DE2, levels = c("significant down", "not significant", "significant up"))
res_W2Relapse_vs_W2Cure_df$DE2_labels <- ifelse(res_W2Relapse_vs_W2Cure_df$DE2 != "not significant", res_W2Relapse_vs_W2Cure_df$gene, NA)

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

my_volcano <- res_W2Relapse_vs_W2Cure_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE2, label = DE2_labels, text = gene)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "DESeq2 W2 relapse vs W2 cure Log2Fold=2 >60%TxnCov") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- res_W2Relapse_vs_W2Cure_df %>% filter(DE2 == "significant up") %>% nrow()
text_down <- res_W2Relapse_vs_W2Cure_df %>% filter(DE2 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano



################################################
############## METARESULTS DESeq ###############

# Do a p-value correction
`DESeq_W2.relapse.ComparedTo.W2.cure` <- read.delim("Data/DE_Run1to3_60TxnCov/W2.cure_vs_W2.relapse/W2_relapse.MTb.DESeq.Ratio.txt") %>% 
  filter(str_detect(GENE_ID, "^Rv\\d+.*"))
DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE <- p.adjust(DESeq_W2.relapse.ComparedTo.W2.cure$PVALUE, method = "fdr")

# Columns for Log2Fold > 2
DESeq_W2.relapse.ComparedTo.W2.cure$DE2 <- ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$LOG2FOLD < -2 & DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE < 0.05, "significant down", ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$LOG2FOLD > 2 & DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE < 0.05, "significant up", "not significant"))
DESeq_W2.relapse.ComparedTo.W2.cure$DE2 <- factor(DESeq_W2.relapse.ComparedTo.W2.cure$DE2, levels = c("significant down", "not significant", "significant up"))
DESeq_W2.relapse.ComparedTo.W2.cure$DE2_labels <- ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$DE2 != "not significant", DESeq_W2.relapse.ComparedTo.W2.cure$GENE_NAME, NA)

my_volcano <- DESeq_W2.relapse.ComparedTo.W2.cure %>%
  ggplot(aes(x = LOG2FOLD, y = -log10(FDR_PVALUE), col = DE2, label = DE2_labels, text = GENE_NAME)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "METARESULTS DESeq2 W2 relapse vs W2 cure Log2Fold=2 >60%TxnCov") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- DESeq_W2.relapse.ComparedTo.W2.cure %>% filter(DE2 == "significant up") %>% nrow()
text_down <- DESeq_W2.relapse.ComparedTo.W2.cure %>% filter(DE2 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano


################################################
######## COMPARE METARESULTS AND DESEQ2 ########

# res_W2Relapse_vs_W2Cure_df vs DESeq_W2.relapse.ComparedTo.W2.cure

setdiff(res_W2Relapse_vs_W2Cure_df$DE2_labels, DESeq_W2.relapse.ComparedTo.W2.cure$DE2_labels)
write.csv(res_W2Relapse_vs_W2Cure_df, "Data/DESeq2_GSEA_Output/res_W2Relapse_vs_W2Cure_df5.csv")
write.csv(DESeq_W2.relapse.ComparedTo.W2.cure, "Data/DESeq2_GSEA_Output/DESeq_W2.relapse.ComparedTo.W2.cure.csv")

################################################
############### METARESULTS ALL ################

my_data <- W2.relapse.ComparedTo.W2.cure

GoodSputum60_RawReadsf <- GoodSamples60_RawReadsf %>% column_to_rownames("X") %>%
  dplyr::select(contains("W"))

# Keep genes that have at least 10 counts in at least 25 samples... Not sure what the best cutoff here is...
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering: Says to use the minimum number of samples in a group
GenesToKeepList <- rownames(GoodSputum60_RawReadsf)[rowSums(GoodSputum60_RawReadsf >= 10) >= 52]

my_data_f <- my_data %>% filter(GENE_ID %in% GenesToKeepList)

# Do a p-value correction
my_data_f$FDR_PVALUE <- p.adjust(my_data_f$AVG_PVALUE, method = "fdr")

# Columns for Log2Fold > 2
DESeq_W2.relapse.ComparedTo.W2.cure$DE2 <- ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$LOG2FOLD < -2 & DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE < 0.05, "significant down", ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$LOG2FOLD > 2 & DESeq_W2.relapse.ComparedTo.W2.cure$FDR_PVALUE < 0.05, "significant up", "not significant"))
DESeq_W2.relapse.ComparedTo.W2.cure$DE2 <- factor(DESeq_W2.relapse.ComparedTo.W2.cure$DE2, levels = c("significant down", "not significant", "significant up"))
DESeq_W2.relapse.ComparedTo.W2.cure$DE2_labels <- ifelse(DESeq_W2.relapse.ComparedTo.W2.cure$DE2 != "not significant", DESeq_W2.relapse.ComparedTo.W2.cure$GENE_NAME, NA)


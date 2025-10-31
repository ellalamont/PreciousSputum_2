
# Can I get all the way through the way people normally do?

source("Import_data.R") # GoodSamples80_RawReadsf, GoodSamples80_pipeSummary


################################################
################ LOAD PACKAGES #################

# https://github.com/JulioLeonIncio/Tutorial-from-DEseq2-to-GSEA-in-R/blob/main/DEseq2_to_GSEA_JL.Rmd
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("DESeq2","sva","fgsea","clusterProfiler","GSEABase","tidyverse","pheatmap")
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


################################################
################ ORGANIZE DATA #################

# Start with GoodSamples80_RawReadsf
my_RawCounts <- GoodSamples80_RawReadsf %>% # rename("GENE_ID" = "X")
  column_to_rownames("X") %>%
  mutate_if(is.numeric, ~ as.integer(round(., 0)))

# Remove any gene with less than 10 reads across all samples (require that for every gene: at least 1 of 6 samples must have counts greater than 10) https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/
# which(rowSums(my_RawCounts >= 10) >= 1)
# my_RawCounts <- my_RawCounts[which(rowSums(my_RawCounts >= 10) >= 1),] # Now only 4025 genes


# Start with GoodSamples80_pipeSummary
my_metadata <- GoodSamples80_pipeSummary %>% 
  mutate(batch = case_when(
    Run == "ProbeTest5" ~ 1,
    Run == "PredictTB_Run1" ~ 2,
    Run == "PredictTB_Run2" ~ 3,
    Run == "PredictTB_Run2.5" ~ 4,
    Run == "PredictTB_Run3" ~ 5)) %>%
  dplyr::select(SampleID, SampleID2, Type2, batch) %>%
  rename("sample" = "SampleID2", "condition" = "Type2") %>%
  column_to_rownames("sample")
my_metadata$condition <- gsub(" ", "_", my_metadata$condition)
my_metadata$condition <- gsub("\\(|\\)", "", my_metadata$condition)
my_metadata$condition <- factor(my_metadata$condition)
my_metadata$batch <- factor(my_metadata$batch)

# Ensure sample names line up
stopifnot(all(colnames(my_RawCounts) %in% rownames(my_metadata)))
my_metadata <- my_metadata[colnames(my_RawCounts), ]
stopifnot(all(colnames(my_RawCounts) == rownames(my_metadata)))

# PCA to visualize data before correction
summary(colSums(my_RawCounts))                   # total reads per sample
hist(log10(colSums(my_RawCounts)+1))             # distribution of library sizes

################################################
################## RUN DESEQ2 ##################

####### TRY WITHOUT COMBAT SEQ ######
# Build DESeq
dds <- DESeqDataSetFromMatrix(countData = my_RawCounts,
                                     colData = my_metadata,
                                     design = ~ condition)
dds <- DESeq(dds)
res_W0Cure_vs_W2Cure <- results(dds, contrast = c("condition", "W2_sputum_cure", "W0_sputum_cure"))

res_W0Cure_vs_W2Cure_Ordered <- res_W0Cure_vs_W2Cure[order(res_W0Cure_vs_W2Cure$pvalue),]
summary(res_W0Cure_vs_W2Cure)

# Save output
res_W0Cure_vs_W2Cure_df <- as.data.frame(res_W0Cure_vs_W2Cure) %>% rownames_to_column("gene")
# res_W0Cure_vs_W2Cure_df$FDR_PVALUE <- p.adjust(res_W0Cure_vs_W2Cure_df$pvalue, method = "fdr")

# Columns for Log2Fold > 2
res_W0Cure_vs_W2Cure_df$DE2 <- ifelse(res_W0Cure_vs_W2Cure_df$log2FoldChange < -2 & res_W0Cure_vs_W2Cure_df$pvalue < 0.05, "significant down", ifelse(res_W0Cure_vs_W2Cure_df$log2FoldChange > 2 & res_W0Cure_vs_W2Cure_df$pvalue < 0.05, "significant up", "not significant"))
res_W0Cure_vs_W2Cure_df$DE2 <- factor(res_W0Cure_vs_W2Cure_df$DE2, levels = c("significant down", "not significant", "significant up"))
res_W0Cure_vs_W2Cure_df$DE2_labels <- ifelse(res_W0Cure_vs_W2Cure_df$DE2 != "not significant", res_W0Cure_vs_W2Cure_df$gene, NA)

# Volcano plot
# ggplot(res_W0Cure_vs_W2Cure_df, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(color = padj < 0.05)) +
#   theme_bw() + scale_color_manual(values = c("grey", "red")) +
#   labs(title = "W0_cure vs Broth")

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

my_volcano <- res_W0Cure_vs_W2Cure_df %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = DE2, label = DE2_labels, text = gene)) + # text is for plotly, could be GENE_ID
  geom_point(alpha = 0.7) + 
  labs(title = "DESeq2 W0_cure vs Broth Log2Fold=2") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel # Changed from 3 to 4
  
  # Need it this way so the colors aren't messed up by not having significant up or down
  # scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) + 
  scale_color_manual(values = c(`significant down` = "#00AFBB", `not significant` = "grey", `significant up` = "#bb0c00"))
plot_build <- ggplot_build(my_volcano)
y_max <- max(plot_build$layout$panel_scales_y[[1]]$range$range)
x_max <- max(plot_build$layout$panel_scales_x[[1]]$range$range)
x_min <- min(plot_build$layout$panel_scales_x[[1]]$range$range)
text_up <- res_W0Cure_vs_W2Cure_df %>% filter(DE2 == "significant up") %>% nrow()
text_down <- res_W0Cure_vs_W2Cure_df %>% filter(DE2 == "significant down") %>% nrow()
my_volcano_annotated <- my_volcano +
  annotate("label", x = (x_max+1)/2, y = y_max - 0.1, label = paste0(text_up, " genes"), color = "#bb0c00", fontface = "bold", fill = "transparent", label.size = 0.3) + 
  annotate("label", x = (x_min-1)/2, y = y_max - 0.1, label = paste0(text_down, " genes"), color = "#00AFBB", fontface = "bold", fill = "transparent", label.size = 0.3)
final_volcano <- my_volcano_annotated + my_plot_themes
final_volcano

################################################
##################### GSEA #####################

# Remove NAs 
res <- res_W0Cure_vs_W2Cure_df %>%
  filter(!is.na(padj)) %>% 
  arrange(desc(log2FoldChange))

# Create a ranked vector (using stat but think could also use Log2FoldChange)
ranks <- res$log2FoldChange
names(ranks) <- res$gene
ranks <- sort(ranks, decreasing = TRUE)

# Optional: filter out NAs or duplicated genes
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[!duplicated(names(ranks))]

# --- 6. Load your custom gene sets ---
# CSV format: Gene,GeneSet
# gene_sets <- read.csv("Data/GeneSet_Data/TAR_Poonawala2024_GeneSets.csv")
# # Convert to list format for enrichment testing
# custom_list <- gene_sets %>%
#   group_by(GeneSet) %>%
#   summarize(genes = list(unique(Gene))) %>%
#   deframe()
custom_list <- allGeneSetList[["EllaGeneSets_2025.10.24"]]

# --- 8. Run fgsea with your gene sets ---
fgsea_res <- fgsea(
  pathways = custom_list,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

fgsea_res_tidy <- fgsea_res %>%
  arrange(padj) %>%
  as_tibble()

# 10. Visualize top pathways ---
topPathwaysUp <- fgsea_res_tidy %>%
  filter(ES > 0) %>%
  arrange(padj) %>%
  head(10)

topPathwaysDown <- fgsea_res_tidy %>%
  filter(ES < 0) %>%
  arrange(padj) %>%
  head(10)

ggplot(topPathwaysUp, aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(title = "Top Upregulated Pathways", x = "Pathway", y = "Normalized Enrichment Score")

ggplot(topPathwaysDown, aes(reorder(pathway, NES), NES)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Downregulated Pathways", x = "Pathway", y = "Normalized Enrichment Score")












### OLD BELOW ###

# dds_raw <- DESeqDataSetFromMatrix(countData = my_RawCounts,
#                                   colData = my_metadata,
#                                   design = ~ condition)
# dds_raw


# combat_counts <- ComBat_seq(counts = as.matrix(my_RawCounts),
#                             batch = my_metadata$batch)
# 
# # Build DESeq on ComBat-corrected counts
# dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
#                                      colData = my_metadata,
#                                      design = ~ condition)   # now omit batch because already corrected
# dds_combat <- DESeq(dds_combat)
# res_W0Cure_vs_Ra <- results(dds_combat, contrast = c("condition", "W0 sputum (cure)", "Broth"))
# 
# res_W0Cure_vs_Ra_Ordered <- res_W0Cure_vs_Ra[order(res_W0Cure_vs_Ra$padj),]
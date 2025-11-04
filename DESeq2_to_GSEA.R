# DESeq2 to GSEA with Raw Reads and custom gene sets
# E. Lamont
# 11/3/25


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
# 11/4/25: Changed to be >60% txn coverage!

# Start with GoodSamples60_RawReadsf
my_RawCounts <- GoodSamples60_RawReadsf %>% # rename("GENE_ID" = "X")
  column_to_rownames("X") %>%
  mutate_if(is.numeric, ~ as.integer(round(., 0)))

# Start with GoodSamples60_pipeSummary
my_metadata <- GoodSamples60_pipeSummary %>% 
  mutate(batch = case_when(
    Run == "ProbeTest5" ~ 1,
    Run == "PredictTB_Run1" ~ 2,
    Run == "PredictTB_Run2" ~ 3,
    Run == "PredictTB_Run2.5" ~ 4,
    Run == "PredictTB_Run3" ~ 5)) %>% # This is for batch correction which I haven't done
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

################################################
################## RUN DESEQ2 ##################
# Without CombatSeq

# Build DESeq
dds <- DESeqDataSetFromMatrix(countData = my_RawCounts,
                                     colData = my_metadata,
                                     design = ~ condition)
dds <- DESeq(dds)
res_W0Cure_vs_W2Cure <- results(dds, contrast = c("condition", "W2_sputum_cure", "W0_sputum_cure"))

summary(res_W0Cure_vs_W2Cure)

# Save output
res_W0Cure_vs_W2Cure_df <- as.data.frame(res_W0Cure_vs_W2Cure) %>% rownames_to_column("gene")

# Columns for Log2Fold > 2
res_W0Cure_vs_W2Cure_df$DE2 <- ifelse(res_W0Cure_vs_W2Cure_df$log2FoldChange < -2 & res_W0Cure_vs_W2Cure_df$padj < 0.05, "significant down", ifelse(res_W0Cure_vs_W2Cure_df$log2FoldChange > 2 & res_W0Cure_vs_W2Cure_df$padj < 0.05, "significant up", "not significant"))
res_W0Cure_vs_W2Cure_df$DE2 <- factor(res_W0Cure_vs_W2Cure_df$DE2, levels = c("significant down", "not significant", "significant up"))
res_W0Cure_vs_W2Cure_df$DE2_labels <- ifelse(res_W0Cure_vs_W2Cure_df$DE2 != "not significant", res_W0Cure_vs_W2Cure_df$gene, NA)


################################################
################# VOLCANO PLOT #################

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
  labs(title = "DESeq2 W2_cure vs W0_cure Log2Fold=2") + 
  geom_vline(xintercept = c(-2,2), col = "grey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), col = "grey", linetype = "dashed") + 
  geom_text_repel(max.overlaps = 10, size = 4) +  # Can do geom_text_repel or geom_label_rebel 
  # Need it this way so the colors aren't messed up by not having significant up or down
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

# Create a ranked vector (using stat or Log2FoldChange)
ranks <- res$log2FoldChange
names(ranks) <- res$gene
ranks <- sort(ranks, decreasing = TRUE)

# Filter out NAs or duplicated genes
ranks <- ranks[!is.na(ranks)]
ranks <- ranks[!duplicated(names(ranks))]

# Load custom gene sets. Needs to have columns Gene, GeneSet. Will take from allGeneSetList
custom_list <- allGeneSetList[["EllaGeneSets_2025.10.24"]]
EllaGeneSets_2025.10.24 <- read.csv("Data/GeneSet_Data/EllaGeneSets_2025.10.24.csv") # Also need the csv to get the grouping for plotting with facets

# Run fgsea
fgsea_res <- fgsea(
  pathways = custom_list,
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

################################################
################### DOT PLOT ###################
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
  left_join(EllaGeneSets_2025.10.24 %>% # Add the Group names
              rename(pathway = GeneSet) %>%
              dplyr::select(pathway, Group), by = "pathway")

# Make the bubble plot
my_bubblePlot <- fgsea_res2 %>%
  mutate(Group_wrapped = str_wrap(Group, width = 19)) %>%
  mutate(Group_wrapped = case_when(Group_wrapped == "Ribosomal proteins" ~ "Ribosomal\nproteins", Group_wrapped == "Hypoxia related" ~ "Hypoxia\nrelated", TRUE ~ Group_wrapped)) %>%
  filter(Group != "Toxin/Antitoxin") %>% # Remove this because I don't think its interesting
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
#        file = paste0("W2CureVsW0Cure_EllaGeneSets_2024.10.24_DESeq2_GSEA", ".pdf"),
#        path = "Figures/Bubbles/DESeq2_GSEA/EllaGeneSets",
#        width = 8, height = 8, units = "in")


################################################
######## RUN DESEQ2 - W2 RELAPSE vs CURE #######
# Without CombatSeq

# Build DESeq
dds <- DESeqDataSetFromMatrix(countData = my_RawCounts,
                              colData = my_metadata,
                              design = ~ condition)
dds <- DESeq(dds)
res_W2Relapse_vs_W2Cure <- results(dds, contrast = c("condition", "W2_sputum_relapse", "W2_sputum_cure"))


summary(res_W2Relapse_vs_W2Cure)

# Save output
res_W2Relapse_vs_W2Cure_df <- as.data.frame(res_W2Relapse_vs_W2Cure) %>% rownames_to_column("gene")
# write.csv(res_W2Relapse_vs_W2Cure_df, "Data/DESeq2_GSEA_Output/res_W2Relapse_vs_W2Cure_df.csv")

# Columns for Log2Fold > 2
res_W2Relapse_vs_W2Cure_df$DE2 <- ifelse(res_W2Relapse_vs_W2Cure_df$log2FoldChange < -2 & res_W2Relapse_vs_W2Cure_df$padj < 0.05, "significant down", ifelse(res_W2Relapse_vs_W2Cure_df$log2FoldChange > 2 & res_W2Relapse_vs_W2Cure_df$padj < 0.05, "significant up", "not significant"))
res_W2Relapse_vs_W2Cure_df$DE2 <- factor(res_W2Relapse_vs_W2Cure_df$DE2, levels = c("significant down", "not significant", "significant up"))
res_W2Relapse_vs_W2Cure_df$DE2_labels <- ifelse(res_W2Relapse_vs_W2Cure_df$DE2 != "not significant", res_W2Relapse_vs_W2Cure_df$gene, NA)


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
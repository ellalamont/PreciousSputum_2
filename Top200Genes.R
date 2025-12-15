# Look at the top 200 genes
# E. Lamont
# 12/15/25

source("Import_data.R") 

# Tried with Txn Coverage >= 20% and nothing really showed up... so go back to 60% threshold

###########################################################
##################### COLLECT DATA ########################

# Since we are using the top 200 genes

SputumSamples60_pipeSummary <- my_pipeSummary %>% filter(SampleID2 %in% SputumSampleList60) %>%
  mutate(Outcome2 = case_when(Outcome == "Prob Relapse" ~ "Relapse", TRUE ~ Outcome))

###########################################################
############## W0: GENES ENRICHED IN RELAPSE ##############

# Get a clean metadata for the W0 samples
SputumSamples60_metadata <- SputumSamples60_pipeSummary %>% 
  filter(Type == "Week 0 sputum") %>%
  select(SampleID2, Outcome2)

# Get a clean tpm for the W0 samples
W0SputumSamples60_tpmf <- All_tpm_f %>% dplyr::select(all_of(SputumSamples60_metadata %>% pull(SampleID2)))

# Check everything is good
stopifnot(all(SputumSamples60_metadata$SampleID2 == colnames(W0SputumSamples60_tpmf)))

# Get the top 200 genes per sample
GetTopGenes_Function <- function(tpm_df, n=200) {
  names(sort(tpm_df, decreasing=T))[1:n]
}
top200Genes_bySample <- apply(W0SputumSamples60_tpmf, 2, GetTopGenes_Function)

# Split samples by outcome
cure_samples <- SputumSamples60_metadata %>% filter(Outcome2 == "Cure") %>% pull(SampleID2)
relapse_samples <- SputumSamples60_metadata %>% filter(Outcome2 == "Relapse") %>% pull(SampleID2)

# Count how many times each gene appears in the cure or relapse top 200 gene lists
countGenes_Function <- function(samples) {
  unlist(top200Genes_bySample[, samples]) %>%
    base::table() %>%
    as.data.frame() %>%
    rename(Gene = ".", Count = Freq)
}
cure_counts <- countGenes_Function(cure_samples)
relapse_counts <- countGenes_Function(relapse_samples)

# Convert counts into percent of samples
# Currently want the genes that are in >=70% of the samples
cure_counts <- cure_counts %>%
  mutate(CurePercent = Count/length(cure_samples) * 100)
relapse_counts <- relapse_counts %>%
  mutate(RelapsePercent = Count/length(relapse_samples) * 100)

# Combine and calculate differences??
topGeneComparison <- full_join(
  cure_counts %>% select(Gene, CurePercent), 
  relapse_counts %>% select(Gene, RelapsePercent), by = "Gene") %>%
  replace_na(list(CurePercent = 0, RelapsePercent = 0)) %>%
  mutate(Difference = RelapsePercent - CurePercent)

# Define the genes enriched in relapse samples
# Top genes in >=70% of relapse samples and enriched vs cure samples
relapse_enriched <- topGeneComparison %>%
  filter(RelapsePercent >= 70 & CurePercent <= 50) %>%
  arrange(desc(RelapsePercent), desc(Difference))
# genes that are in 70% of relapse and in less than 50% of cure cases:
# Rv0684  Rv0145  Rv3458c Rv1375  Rv2660c Rv3557c


# # Visualization - Pretty weird right now
# presence_matrix <- sapply(colnames(W0SputumSamples60_tpmf),
#                           function(s) rownames(W0SputumSamples60_tpmf) %in% top200Genes_bySample[,s])
# presence_matrix_num <- apply(presence_matrix, 2, as.numeric) # 1=Gene in top 200
# rownames(presence_matrix_num) <- rownames(W0SputumSamples60_tpmf)
# My_annotation_col <- SputumSamples60_metadata %>% column_to_rownames("SampleID2") # Need to do this so the pheatmap works
# pheatmap(presence_matrix_num[relapse_enriched$Gene,],
#          annotation_col = My_annotation_col,
#          show_rownames = T,
#          cluster_rows = F)

###########################################################
############## W2: GENES ENRICHED IN RELAPSE ##############

# Get a clean metadata for the W0 samples
SputumSamples60_metadata <- SputumSamples60_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>% 
  filter(Outcome2 != "Failure") %>%
  select(SampleID2, Outcome2)

# Get a clean tpm for the W0 samples
W2SputumSamples60_tpmf <- All_tpm_f %>% dplyr::select(all_of(SputumSamples60_metadata %>% pull(SampleID2)))

# Check everything is good
stopifnot(all(SputumSamples60_metadata$SampleID2 == colnames(W2SputumSamples60_tpmf)))

# Get the top 200 genes per sample
GetTopGenes_Function <- function(tpm_df, n=200) {
  names(sort(tpm_df, decreasing=T))[1:n]
}
top200Genes_bySample <- apply(W2SputumSamples60_tpmf, 2, GetTopGenes_Function)

# Split samples by outcome
cure_samples <- SputumSamples60_metadata %>% filter(Outcome2 == "Cure") %>% pull(SampleID2)
relapse_samples <- SputumSamples60_metadata %>% filter(Outcome2 == "Relapse") %>% pull(SampleID2)

# Count how many times each gene appears in the cure or relapse top 200 gene lists
countGenes_Function <- function(samples) {
  unlist(top200Genes_bySample[, samples]) %>%
    base::table() %>%
    as.data.frame() %>%
    rename(Gene = ".", Count = Freq)
}
cure_counts <- countGenes_Function(cure_samples)
relapse_counts <- countGenes_Function(relapse_samples)

# Convert counts into percent of samples
# Currently want the genes that are in >=70% of the samples
cure_counts <- cure_counts %>%
  mutate(CurePercent = Count/length(cure_samples) * 100)
relapse_counts <- relapse_counts %>%
  mutate(RelapsePercent = Count/length(relapse_samples) * 100)

# Combine and calculate differences??
topGeneComparison <- full_join(
  cure_counts %>% select(Gene, CurePercent), 
  relapse_counts %>% select(Gene, RelapsePercent), by = "Gene") %>%
  replace_na(list(CurePercent = 0, RelapsePercent = 0)) %>%
  mutate(Difference = RelapsePercent - CurePercent)

# Define the genes enriched in relapse samples
# Top genes in >=70% of relapse samples and enriched vs cure samples
relapse_enriched <- topGeneComparison %>%
  filter(RelapsePercent >= 70 & CurePercent <= 50) %>%
  arrange(desc(RelapsePercent), desc(Difference))
# genes that are in 70% of relapse and in less than 50% of cure cases:
# Rv0581  Rv1884c Rv2021c Rv1080c
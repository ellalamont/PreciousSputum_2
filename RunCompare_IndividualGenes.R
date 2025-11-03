# Compare individual genes between samples 
# After meeting with David 11/3/25
# E. Lamont

# Are the genes with <10 TPM the same across all samples? Between identical samples, W2 samples and between W0 and W2 samples


source("Import_data.R") 


###########################################################
########### W2 HIGHLY/LOWLY EXPRESSED GENES ###############

# I guess look at the ones above 80% and above 60%?

W2_80txnCov_tpmf <- GoodSamples80_tpmf %>% select(contains("W2"))

# I want a list of genes that are <10 TPM in each sample.. maybe better to pull from csv and put into interativenn

write.csv(W2_80txnCov_tpmf, "Data/Interactivenn/W2_80txnCov_tpmf.csv")

# After looking at interactivenn, not much overlap for the lowly expressed samples, but 93/200 top 200 TPMf samples are the same, see what these are

# Grab the top 200 genes for all of these samples and make a list of lists
Top200Genes <- lapply(names(W2_80txnCov_tpmf), function(sample) {
  W2_80txnCov_tpmf %>%
    arrange(desc(.data[[sample]])) %>%
    slice_head(n = 200) %>%
    rownames()})
names(Top200Genes) <- names(W2_80txnCov_tpmf)

# Find the genes that are in the top 200 genes in all the samples
W2_80TxnCov_TopGenes_Intersection <- Reduce(intersect, Top200Genes) # 93, which is what interactivenn also gave



###########################################################
########### W0 HIGHLY/LOWLY EXPRESSED GENES ###############

# Won't be able to do this in interactivenn because there are so many more samples....

W0_80txnCov_tpmf <- GoodSamples80_tpmf %>% select(contains("W0"))

# Grab the top 200 genes for all of these samples and make a list of lists
Top200Genes <- lapply(names(W0_80txnCov_tpmf), function(sample) {
  W0_80txnCov_tpmf %>%
    arrange(desc(.data[[sample]])) %>%
    slice_head(n = 200) %>%
    rownames()})
names(Top200Genes) <- names(W0_80txnCov_tpmf)

# Find the genes that are in the top 200 genes in all the samples
W0_80TxnCov_TopGenes_Intersection <- Reduce(intersect, Top200Genes) # Only 18 here

# NEXT: Are there 18 genes in the W2 list?



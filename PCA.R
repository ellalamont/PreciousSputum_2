# PCA on all the biological samples passing filter
# E. Lamont
# 9/26/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R") # To get GoodSamples_tpmf and GoodSamples_pipeSummary 


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )

# Labelled Colors
my_fav_colors <- c(`W0 sputum (cure)` = "#0072B2", `W0 sputum (relapse)` = "green4", `Broth`= "#999999", `W2 sputum (cure)` = "#6A3D9A", `W2 sputum (relapse)` = "#E69F00", `THP1 spiked` = "red")
my_fav_colors2 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00") 
# Labelled Shapes
my_fav_shapes <- c(`Sputum` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)
my_fav_shapes2 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25)

###########################################################
################## PCA BIOL with BROTH ####################
# Passing filter is >1,000,000 genomic reads and >80% genes with at least 10 reads, already subsetted in Import_data.R

# Convert gene column to rownames
my_tpm <- GoodBiolSamples_tpm %>% column_to_rownames(var = "X")

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 34.1% of variance
summary_PCA[2,1] # PC2 explains 12.4% of variance
summary_PCA[3,1] # PC3 explains 9.1% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = Lineage), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors2) +  
  scale_shape_manual(values = my_fav_shapes2) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(# title = "PCA: >1M reads and >80% genes with at least 10 reads",
       # subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_GoodSamples_2.pdf"),
       path = "Figures_preNonCodingRemoval/PCA",
       width = 6, height = 6, units = "in")


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Type3# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D


###########################################################
############# WHAT GENES IS EACH PC MADE OF ###############

# View(my_PCA$rotation)

My_PCA_rotation <- as.data.frame(my_PCA$rotation) %>%
  rownames_to_column("Gene")

# Plot the PC rotation values
PC_to_plot <- "PC1"
top_n_genes <- 20

# Get the top Gene contributors for the PCA (positive or negative values)
top_genes_PC1 <- My_PCA_rotation %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:50) %>%
  select(Gene, PC1)
top_genes_PC2 <- My_PCA_rotation %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:50) %>%
  select(Gene, PC2)

# Add information about the genes 
load("Data/MTb.MapSet.rda")
my_geneInfo <- mapSet[["geneMap"]] %>% select(GENE_ID, NAME, PRODUCT)
top_genes_PC1 <- inner_join(top_genes_PC1 %>% rename(GENE_ID = Gene), 
                           my_geneInfo, by = "GENE_ID")
top_genes_PC2 <- inner_join(top_genes_PC2 %>% rename(GENE_ID = Gene), 
                            my_geneInfo, by = "GENE_ID")

# Save these
write.csv(top_genes_PC1, "Figures/PCA/TPM_GoodSamples_TopGenes_PC1.csv")
write.csv(top_genes_PC2, "Figures/PCA/TPM_GoodSamples_TopGenes_PC2.csv")


###########################################################
############# PCA BIOL with BROTH FILTERED ################
# Filtered meaning the non coding genes have been removed

# Convert gene column to rownames
my_tpm <- GoodBiolSamples_tpm_f # %>% column_to_rownames(var = "X")

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 34.7% of variance
summary_PCA[2,1] # PC2 explains 12.6% of variance
summary_PCA[3,1] # PC3 explains 9.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type, shape = Type)) + 
  geom_point(aes(fill = Type, shape = Type), size = 5, alpha = 0.7, stroke = 0.8) +
  geom_text_repel(aes(label = Lineage), size = 2.5, max.overlaps = Inf) + 
  scale_fill_manual(values = my_fav_colors) +  
  scale_shape_manual(values = my_fav_shapes) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM; Rvnc removed",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
# ggsave(PCA_tpm_1,
#        file = paste0("TPMf_GoodSamples_1.pdf"),
#        path = "Figures/PCA",
#        width = 8, height = 5, units = "in")
# ggplotly(PCA_tpm_1)


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Type# , 
                  # colors = c12,
                  # text = ~Replicate
)
PCA_3D



###########################################################
################ PCA BIOL with BROTH and Rv ###############

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027", `IndigoRv` = "black", `NoDrugRv` = "#E7298A") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23, `IndigoRv` = 23, `NoDrugRv` = 23)

# Passing filter is >1,000,000 genomic reads and >80% genes with at least 10 reads, already subsetted in Import_data.R

# Convert gene column to rownames
my_tpm <- GoodBiolSamples_w_Rv_tpm %>% column_to_rownames(var = "X")

# Remove the Ra, and others
# my_tpm <- my_tpm %>% select(-contains("BQ"), -contains("HN878"), -contains("Cav"))

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 28.2% of variance
summary_PCA[2,1] # PC2 explains 15.6% of variance
summary_PCA[3,1] # PC3 explains 10.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
# ggsave(PCA_tpm_1,
#        file = paste0("TPM_wRv_2.pdf"),
#        path = "Figures_preNonCodingRemoval/PCA/wRv",
#        width = 8, height = 5, units = "in")

###########################################################
####################### PCA with Rv #######################
# 9/11/25: Just the pH7 Rv, not Ra

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027", `IndigoRv` = "black", `NoDrugRv` = "#E7298A") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23, `IndigoRv` = 23, `NoDrugRv` = 23)

# Passing filter is >1,000,000 genomic reads and >80% genes with at least 10 reads, already subsetted in Import_data.R

# Convert gene column to rownames
my_tpm <- GoodBiolSamples_w_Rv_tpm %>% column_to_rownames(var = "X")

# Remove the Ra, and others
# my_tpm <- my_tpm %>% select(-contains("BQ"), -contains("HN878"), -contains("Cav"))

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Remove the Ra and other pH Rvs
my_tpm_t2_subset <- my_tpm_t2 %>% filter(!rownames(my_tpm_t2) %in% c("H37Ra_Broth_4_S7", "H37Ra_Broth_5_S8", "H37Ra_Broth_6_S9", "Rv_pH_8.3_R1", "Rv_pH_8.3_R2", "Rv_pH_5.7_R1", "Rv_pH_5.7_R2", "INDIGO.NoDrug.Control.E1"))

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2_subset, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 33.9% of variance
summary_PCA[2,1] # PC2 explains 13.3% of variance
summary_PCA[3,1] # PC3 explains 9.8% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
# ggsave(PCA_tpm_1,
#        file = paste0("TPM_wRv_2.pdf"),
#        path = "Figures_preNonCodingRemoval/PCA/wRv",
#        width = 8, height = 5, units = "in")



# Get the top Gene contributors for the PCA (positive or negative values)
top_genes_PC2 <- as.data.frame(my_PCA$rotation) %>%
  rownames_to_column("Gene") %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:20) %>%
  select(Gene, PC2)


###########################################################
####################### PCA TPM BY HAND ###################
# 9/1/25

my_fav_colors3 <- c(`Sputum L4` = "#004C73", `Sputum L2` = "#99CCE8", `Caseum mimic` = "green4", `Broth`= "#999999", `Marmoset` = "#6A3D9A", `Rabbit` = "#E69F00", `Rv7` = "#FC8D59", `Rv8.3` = "#FEE090", `Rv5.7` = "#D73027", `IndigoRv` = "black", `NoDrugRv` = "#E7298A") 
# Labelled Shapes
my_fav_shapes3 <- c(`Sputum L4` = 21, `Sputum L2` = 21, `Caseum mimic` = 22, `Broth`= 23, `Marmoset` = 24, `Rabbit` = 25, `Rv7` = 23, `Rv8.3` = 23, `Rv5.7` = 23, `IndigoRv` = 23, `NoDrugRv` = 23)

# Passing filter is >1,000,000 genomic reads and >80% genes with at least 10 reads, already subsetted in Import_data.R


# Calculate TPM from Raw Reads
# GoodBiolSamples_wRv_RawReads
my_tpm_byHand <- CalculateTPM(GoodBiolSamples_wRv_RawReads)

# Convert gene column to rownames
# my_tpm <- GoodBiolSamples_w_Rv_tpm %>% column_to_rownames(var = "X")

# Remove the Ra, and others
# my_tpm <- my_tpm %>% select(-contains("BQ"), -contains("HN878"), -contains("Cav"))

# Transform the data
my_tpm_t <- as.data.frame(t(my_tpm_byHand))

# Remove columns that are all zero so the scale works for prcomp
my_tpm_t2 <- my_tpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_tpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 32.0% of variance
summary_PCA[2,1] # PC2 explains 13.8% of variance
summary_PCA[3,1] # PC3 explains 8.5% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, BiolSamples_pipeSummary_2, by = "SampleID", )

PCA_tpm_1 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Type3, shape = Type3)) + 
  geom_point(aes(fill = Type3, shape = Type3), size = 5, alpha = 0.7, stroke = 0.8) +
  # geom_text_repel(aes(label = SampleID), size = 2.5) + 
  scale_fill_manual(values = my_fav_colors3) +  
  scale_shape_manual(values = my_fav_shapes3) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "PCA: >1M reads and >80% genes with at least 10 reads",
       subtitle = "TPM from Raw Reads",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
PCA_tpm_1
ggsave(PCA_tpm_1,
       file = paste0("TPM_wRv_3.pdf"),
       path = "Figures_preNonCodingRemoval/PCA/wRv",
       width = 8, height = 5, units = "in")

# Try a heatmap with pheatmap
# E. Lamont 
# 10/29/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Import_data.R")

# Just want the good sputum and broth


###########################################################
###################### PROCESS DATA #######################


my_tpm <- GoodSamples80_tpmf %>% select(-contains("THP1"))


# Filter so there are fewer genes to deal with right now
# Lets filter for 300 just so it's easier for right now!
my_tpm_filtered300 <- my_tpm %>%
  filter(if_all(where(is.numeric), ~ .x >= 300))

my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "W2 sputum (cure)" = "green4",
            "W2 sputum (relapse)" = "#6A3D9A", 
            "Broth" = "#999999")
)

###########################################################
######################## PHEATMAP #########################

# pheatmap(All_tpm_matrix[1:10,], scale = "row")

testing <- my_tpm_filtered300 %>% subset(rownames(my_tpm_filtered300) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]]) # Guess this doesn't need to be a matrix

# pheatmap(testing, scale = "row")



# Grab the columns needed to give colors
Color_annotation_df <- GoodSamples80_pipeSummary %>%
  filter(SampleID2 %in% colnames(testing)) %>%
  select(SampleID2, Type2) %>%
  column_to_rownames("SampleID2")

# Reorder annotation rows to match columns of tpm file
Color_annotation_df <- Color_annotation_df[colnames(testing), , drop = FALSE]

# Define the colors
my_annotation_colors <- list(
  Type2 = c("W0 sputum (cure)" = "#0072B2",
            "W0 sputum (relapse)" = "red", 
            "Caseum mimic" = "green4",
            "Marmoset" = "#6A3D9A", 
            "Rabbit" = "#E69F00", 
            "Broth" = "#999999")
)

pheatmap(my_tpm_filtered300, 
         # annotation_col = annotation_df, 
         # annotation_colors = my_annotation_colors,
         scale = "row")





###########################################################
######################## ALL DATA #########################

pheatmap(All_tpm2, 
         annotation_col = annotation_df, 
         annotation_colors = my_annotation_colors,
         scale = "row")


All_tpm3 <- All_tpm2[rowSums(All_tpm2 == 0) != ncol(All_tpm2), ]

All_tpm3_matrix <- All_tpm3 %>% 
  # rename("W0_250754" = "S_250754",
  #        "W0_355466" = "S_355466",
  #        "W0_503557" = "S_503557",
  #        "W2_503937" = "S_503937",
  #        "W2_575533" = "S_575533_MtbrRNA",
  #        "W2_577208" = "S_577208") %>%
  as.matrix()


###########################################################
############### ALL DATA WITH CLUSTERING ##################
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

# Grab the columns needed to give colors
Color_annotation_df <- pipeSummary_3 %>%
  filter(SampleID %in% colnames(All_tpm3)) %>%
  select(SampleID, Type2) %>%
  column_to_rownames("SampleID")

# Reorder annotation rows to match columns of tpm file
Color_annotation_df <- Color_annotation_df[colnames(All_tpm3), , drop = FALSE]



heatmap_all <- pheatmap(All_tpm3_matrix, 
         annotation_col = Color_annotation_df, 
         annotation_colors = my_annotation_colors,
         scale = "row")
ggsave(heatmap_all,
       file = "heatmap_all_1.pdf",
       path = "Figures/Heatmaps",
       width = 20, height = 20, units = "in")



######## NOTHING BELOW HAS BEEN DONE #########

# Try to get good row annotations based on MTb functional group
# Start with MTb.TB.Phenotypes.AllGeneSets
# Convert to dataframe
Gene_Category <- do.call(rbind, lapply(names(Walter2015GeneSets), function(category) {
  data.frame(Gene = Walter2015GeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")


pheatmap(my_tpm_2_matrix, 
         # annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 8,
         cutree_cols = 5)





# Pull out what the clusters are: 
# use silent = TRUE to suppress the plot
my_heatmap <- pheatmap(my_tpm_2_matrix, 
                       annotation_row = Gene_Category, 
                       annotation_col = my_pipeSummary["Week"],
                       annotation_colors = my_annotation_colors,
                       scale = "row",
                       cutree_rows = 7,
                       cutree_cols = 5,
                       silent = TRUE)
# Extract the row clustering information
row_clusters <- cutree(my_heatmap$tree_row, k = 7)
# Convert to a data frame for easier handling
row_cluster_df <- data.frame(Gene = names(row_clusters), Cluster = row_clusters)
# View the first few rows
head(row_cluster_df)


###########################################################
############ W0 AND BROTH WITH CLUSTERING #################

my_tpm_3_matrix <- my_tpm_2 %>% select(-c("S_503937", "S_577208", "S_575533_MtbrRNA")) %>%
  as.matrix()

Gene_Category <- do.call(rbind, lapply(names(MTb.TB.Phenotypes.AllGeneSets), function(category) {
  data.frame(Gene = MTb.TB.Phenotypes.AllGeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  # filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")

pheatmap(my_tpm_3_matrix, 
         annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 6,
         cutree_cols = 2)











###########################################################
################### TESTING FOR SHINY #####################


allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["microaerophilic: top 25 genes"]]
my_data <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]])
p <- pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_colors = my_annotation_colors,
              scale = "row")
p
heatmap(as.matrix(my_data))


selected_genes <- c("Rv0081", "Rv0494", "Rv2011c", "Rv1473A")
my_data <- my_tpm[rownames(my_tpm) %in% selected_genes, , drop = FALSE]
pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_row = gene_annot["Product"],  # Conditional annotation
              annotation_colors = my_annotation_colors,
              scale = "row", 
              fontsize = 18)




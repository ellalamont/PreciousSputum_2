# Import Bob's differential gene expression data and the metagenesets data
# 10/10/25
# E. Lamont

source("Import_data.R")

# Data is coming from the Lenovo PredictTB_Run2


###########################################################
################### IMPORT BOB's DE DATA ##################

`W0.cure.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W0.cure_vs_Ra/W0_cure.MTb.Meta.JOINED.txt")
`W0.relapse.ComparedTo.W0.cure` <- read.delim("Data/Differential_Expression/W0.cure_vs_W0.relapse/W0_relapse.MTb.Meta.JOINED.txt")
`W2.cure.ComparedTo.W0.cure` <- read.delim("Data/Differential_Expression/W0.cure_vs_W2.cure/W2_cure.MTb.Meta.JOINED.txt")
`W0.relapse.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W0.relapse_vs_Ra/W0_relapse.MTb.Meta.JOINED.txt")
`W2.cure.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W2.cure_vs_Ra/W2_cure.MTb.Meta.JOINED.txt")


###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`W0.cure.ComparedTo.Ra`,
                 `W0.relapse.ComparedTo.W0.cure`, 
                 `W2.cure.ComparedTo.W0.cure`, 
                 `W0.relapse.ComparedTo.Ra`,
                 `W2.cure.ComparedTo.Ra`)

# Make a list of the names
df_names <- c("W0.cure.ComparedTo.Ra",
              "W0.relapse.ComparedTo.W0.cure", 
              "W2.cure.ComparedTo.W0.cure", 
              "W0.relapse.ComparedTo.Ra",
              "W2.cure.ComparedTo.Ra")

# Give the df list the correct df names
names(list_dfs) <- df_names



###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# 10/10/25: Adjusting this so it is all FDR p-values and I have Log2Fold >1 and >2

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Calculate the FDR p-value
  current_df$FDR_PVALUE <- p.adjust(current_df$AVG_PVALUE, method = "fdr")
  
  # Columns for Log2Fold > 1
  current_df$DE1 <- ifelse(current_df$LOG2FOLD < -1 & current_df$FDR_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 1 & current_df$FDR_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE1 <- factor(current_df$DE1, levels = ordered_DE)
  current_df$DE1_labels <- ifelse(current_df$DE1 != "not significant", current_df$GENE_NAME, NA)
  
  # Columns for Log2Fold > 2
  current_df$DE2 <- ifelse(current_df$LOG2FOLD < -2 & current_df$FDR_PVALUE < 0.05, "significant down", ifelse(current_df$LOG2FOLD > 2 & current_df$FDR_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE2 <- factor(current_df$DE2, levels = ordered_DE)
  current_df$DE2_labels <- ifelse(current_df$DE2 != "not significant", current_df$GENE_NAME, NA)

  list_dfs_2[[current_df_name]] <- current_df
}

###########################################################
################# REMOVE NON-CODING GENES #################
# 8/15/25: After talking to DRS, decided to remove all the non-coding genes and all the MT genes, and leave just the coding genes starting with Rv. So need to remove these at the raw read level and calculate new TPM
# The Pathcap people also had issues with ncRNAs: https://www.nature.com/articles/s41598-019-55633-6#Sec8

list_dfs_f <- lapply(list_dfs_2, function(df) {
  df %>% filter(str_detect(GENE_ID, "^Rv\\d+.*"))
})




###########################################################
############# IMPORT BOB's METAGENESETS DATA ##############

`MetaGeneSets_W0.cure.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W0.cure_vs_Ra/W0_cure.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W0.relapse.ComparedTo.W0.cure` <- read.delim("Data/Differential_Expression/W0.cure_vs_W0.relapse/W0_relapse.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2.cure.ComparedTo.W0.cure` <- read.delim("Data/Differential_Expression/W0.cure_vs_W2.cure/W2_cure.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W0.relapse.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W0.relapse_vs_Ra/W0_relapse.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2.cure.ComparedTo.Ra` <- read.delim("Data/Differential_Expression/W2.cure_vs_Ra/W2_cure.MTb.MetaGeneSets.UP.txt")



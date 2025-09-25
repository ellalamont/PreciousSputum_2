source("Import_data.R")


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, size=8, vjust=1, hjust=1),
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
################ COLLECT DATA OF INTEREST #################

# Collecting all the duplicates between the two runs
Run1_Subset <- Run1_tpm %>% select(X, THP1_1e6_1_S67, W0_12043_S32, W0_12082_S45, W0_13094_S46, W0_14136_S52, W0_15072_S48, W0_15083_S50) %>%
  rename(Gene = X, Run1_THP1.Ra1e6 = THP1_1e6_1_S67, Run1_W0_12043 = W0_12043_S32, Run1_W0_12082 = W0_12082_S45,Run1_W0_13094 = W0_13094_S46,  Run1_W0_14136 = W0_14136_S52, Run1_W0_15072 = W0_15072_S48, Run1_W0_15083 = W0_15083_S50)
Run2_Subset <- Run2_tpm %>% select(X, THP1_1e6_1_S44, W0_12043_S45, W0_12082_S46, W0_13094_S47, W0_14136_S50, W0_15072_S48, W0_15083_S49) %>%
  rename(Gene = X, Run2_THP1.Ra1e6 = THP1_1e6_1_S44, Run2_W0_12043 = W0_12043_S45, Run2_W0_12082 = W0_12082_S46, Run2_W0_13094 = W0_13094_S47, Run2_W0_14136 = W0_14136_S50, Run2_W0_15072 = W0_15072_S48, Run2_W0_15083 = W0_15083_S49)

merged_DoubleRun <- merge(Run1_Subset, Run2_Subset, all = T)

# Log10 transform the data
merged_DoubleRun_Log10 <- merged_DoubleRun %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Remove all the non Rv genes and see how much better it gets
merged_DoubleRun_Log10_filtered <- merged_DoubleRun_Log10 %>% 
  filter(grepl("^Rv[0-9]+[A-Za-z]?$", Gene))

# Pivot longer the data
merged_DoubleRun_Log10_t <- merged_DoubleRun_Log10_filtered %>%
  pivot_longer(cols = starts_with("Run"), names_to = "Sample", values_to = "TPM_log10")


# Pull out the sample names 
sample_choices <- merged_DoubleRun_Log10_t %>%
  pull(Sample) %>%
  unique() %>%
  str_remove("^Run[12]_") %>%
  unique()


##################################################################################

# Define UI
ui <- fluidPage(
  titlePanel("Compare Gene Reads Across Runs"),
  sidebarLayout(
    sidebarPanel(
      selectInput("sample_choice", "Choose a Sample:",
                  choices = sample_choices, 
                  selected = sample_choices[1]),
      sliderInput("gene_range", "Gene Window",
                  min = 1, max = length(unique(merged_DoubleRun_Log10_t$Gene)) - 99,
                  value = 1, step = 10, width = "100%")
    ),
    mainPanel(
      plotlyOutput("BarPlot", height = "600px")
    )
  )
)

# Define Server
server <- function(input, output) {
  
  # Choose the two samples and determine gene window
  filtered_data <- reactive({
    df <- merged_DoubleRun_Log10_t %>%
      filter(str_detect(Sample, paste0("Run[12]_", input$sample_choice)))
    unique_genes <- unique(df$Gene)
    gene_start <- as.numeric(input$gene_range)
    gene_end <- min(gene_start + 99, length(unique_genes))
    selected_genes <- unique_genes[gene_start:gene_end]
    df %>% filter(Gene %in% selected_genes)
  })
  

  
  output$BarPlot <- renderPlotly({
    
    p <- filtered_data() %>%
      # head(100) %>% # If you run without this it crashes R
      ggplot(aes(x = Gene, y = TPM_log10, fill = Sample)) +
      geom_bar(position="fill", stat="identity") +
      scale_fill_manual(values = c("#4AAFD5", "#E7A339")) + 
      labs(
        title = paste("Gene hits for sample:", input$sample_choice),
        y = "log10(TPM + 1)",
        x = "Gene"
      ) +
      scale_y_continuous(expand = c(0,0)) + 
      my_plot_themes
    
    ggplotly(p)
  })
  
}

shinyApp(ui, server)
# 11/3/25
# Interactive plot of the corrlation between the average TPM scaled values for THP1 cells spiked with 1e6 cells oF H37Ra.
# Doing this so I can look at where the gene sets are in the correlation!

library(shiny)
library(gmodels)

source("Import_data.R") 

# Plot basics
my_plot_themes <- theme_bw() +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12))

###########################################################
############ TPM_F: COLLECT DATA OF INTEREST ##############

# Genes have been filtered to keep only protein coding Rv genes and then TPM done manually (not Bob's pipeline)

# Log10 transform the data
All_tpmf_Log10 <- All_tpm_f %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Make Gene a column
All_tpmf_Log10 <- All_tpmf_Log10 %>% rownames_to_column("Gene")





# Define UI ----
ui <- fluidPage(
  titlePanel("Correlations to compare runs"),
  
  fluidRow(
    
    column(width = 3,
           textInput("my_GeneID", 
                     label = "Gene ID (Will label gene orange)",
                     value = "Rv..."),
           selectInput("Sample1",
                       label = "Sample 1 (x-axis)",
                       choices = names(All_tpmf_Log10)[!names(All_tpmf_Log10) %in% "Gene"]),
           selectInput("Sample2",
                       label = "Sample 2 (y-axis)",
                       choices = names(All_tpmf_Log10)[!names(All_tpmf_Log10) %in% "Gene"]),
    ),
    
    column(width = 9,
           plotlyOutput("Correlation_plot",
                        width = "90%", height = "600px"),
    ),
  )
  
)


# Define server logic ----
server <- function(input, output) {
  
  # Correlation Plot
  
  output$Correlation_plot <- renderPlotly({
    
    # Add data for labelling a single gene
    single_gene <- All_tpmf_Log10 %>% 
      filter(Gene == input$my_GeneID)
    
    ScatterCorr <- All_tpmf_Log10 %>% 
      ggplot(aes(x = .data[[input$Sample1]], y = .data[[input$Sample2]])) + 
      geom_point(aes(text = Gene), alpha = 0.8, size = 2, color = "black") +
      
      # Add a differently colored point
      geom_point(data = single_gene, color = "orange", aes(label = Gene)) + 
      
      geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
      labs(title = paste0(input$Sample1, " vs ", input$Sample2),
           subtitle = "Pearson correlation",
           x = paste0("Log10(TPM+1) ", input$Sample1), y = paste0("Log10(TPM+1) ", input$Sample2)) + 
      stat_cor(method="pearson") + # add a correlation to the plot
      my_plot_themes
    ScatterCorr
    
    
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
library(shiny)
library(tidyverse)
library(pheatmap)
library(grid)

## 10/14/25 FIXED VERSION ----

# Load your data
source("Import_data.R")

# Plot theme
my_plot_themes <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(angle = 0, size = 12, vjust = 1, hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 20)
  )

# Define annotation colors
my_annotation_colors <- list(
  Week = c("W0" = "#0072B2",  # Blue
           "W2" = "#E66900",  # Orange
           "Broth" = "#999999")  # Grey
)

# Remove THP1 columns
my_tpm <- GoodSamples80_tpmf %>% select(-contains("THP1"))

# Fix week labeling
my_pipeSummary <- GoodSamples80_pipeSummary %>%
  mutate(
    Week = case_when(
      Week == "Week 0" ~ "W0",
      Week == "Week 2" ~ "W2",
      is.na(Week) ~ "Broth",
      TRUE ~ Week
    )
  )

# ---- UI ----
ui <- fluidPage(
  titlePanel("Sputum Pheatmap"),
  
  fluidRow(
    column(
      width = 4,
      
      # Timepoint checkboxes
      checkboxGroupInput(
        "timepoints",
        label = "Select Timepoints",
        choices = c("Broth", "W0", "W2"),
        selected = c("Broth", "W0", "W2")
      ),
      
      # Gene set source
      selectInput(
        "my_GeneSetSource",
        label = "Gene Set Source",
        choices = names(allGeneSetList)
      ),
      
      # Gene set within that source
      selectInput(
        "my_GeneSet",
        label = "Gene Set",
        choices = NULL
      ),
      
      # Manual gene entry
      textInput(
        "manual_genes",
        label = "Enter Rv# (comma separated)",
        placeholder = "Rv1473A, Rv2011c, Rv0494"
      ),
      
      # Clustering options
      numericInput("cutree_rows", "Number of Row Clusters", value = 1, min = 1, step = 1),
      numericInput("cutree_cols", "Number of Column Clusters", value = 1, min = 1, step = 1),
      
      # Scaling
      selectInput("my_scaling", "How to scale", choices = c("row", "column", "none")),
      
      # Display numbers
      checkboxInput("show_numbers", label = "Show Values", value = FALSE),
      
      # Mycobrowser link
      textInput(
        "my_GeneID",
        label = "Pick a gene to link to mycobrowser",
        placeholder = "Rv..."
      ),
      uiOutput("gene_link")
    ),
    
    column(
      width = 8,
      uiOutput("dynamic_pheatmap")
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  # Mycobrowser gene link
  output$gene_link <- renderUI({
    req(input$my_GeneID)
    url <- paste0("https://mycobrowser.epfl.ch/genes/", input$my_GeneID)
    tags$a(
      href = url,
      target = "_blank",
      paste0("View Details of ", input$my_GeneID, " on Mycobrowser")
    )
  })
  
  # Update available gene sets when source changes
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(
      session,
      "my_GeneSet",
      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
      selected = NULL
    )
  })
  
  # Get selected genes
  get_selected_genes <- reactive({
    if (input$manual_genes != "") {
      genes <- unlist(strsplit(input$manual_genes, "[,\\s]+"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
    } else if (!is.null(input$my_GeneSet) &&
               input$my_GeneSet %in% names(allGeneSetList[[input$my_GeneSetSource]])) {
      genes <- allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]]
    } else {
      genes <- character(0)
    }
    return(genes)
  })
  
  # Dynamic plot height
  output$dynamic_pheatmap <- renderUI({
    selected_genes <- get_selected_genes()
    num_genes <- sum(rownames(my_tpm) %in% selected_genes)
    plot_height <- max(400, min(2500, num_genes * 40))
    plotOutput("pheatmap", height = paste0(plot_height, "px"))
  })
  
  # Render the heatmap
  output$pheatmap <- renderPlot({
    
    # Columns for each timepoint
    timepoint_columns <- list(
      "W0" = colnames(my_tpm)[str_detect(colnames(my_tpm), "W0")],
      "W2" = colnames(my_tpm)[str_detect(colnames(my_tpm), "W2")],
      "Broth" = colnames(my_tpm)[str_detect(colnames(my_tpm), "Broth")]
    )
    
    selected_genes <- get_selected_genes()
    
    # Filter rows (genes)
    my_data <- my_tpm[rownames(my_tpm) %in% selected_genes, , drop = FALSE]
    
    # Filter columns (samples)
    if (!is.null(input$timepoints)) {
      selected_cols <- unlist(timepoint_columns[input$timepoints])
      selected_cols <- intersect(selected_cols, colnames(my_data))
      my_data <- my_data[, selected_cols, drop = FALSE]
    }
    
    # Validation checks
    if (ncol(my_data) == 0) {
      showNotification("No columns match the selected timepoints.", type = "error")
      return(NULL)
    }
    
    if (nrow(my_data) < 2) {
      showNotification("At least two valid genes are required for clustering.", type = "error")
      return(NULL)
    }
    
    # Make sure annotation matches column names
    annotation_df <- my_pipeSummary %>%
      select(SampleID, Week) %>%
      filter(SampleID %in% colnames(my_data)) %>%
      column_to_rownames("SampleID")
    
    annotation_df <- annotation_df[colnames(my_data), , drop = FALSE]
    
    # Draw the heatmap
    p <- pheatmap(
      my_data,
      annotation_col = annotation_df,
      annotation_colors = my_annotation_colors,
      scale = input$my_scaling,
      display_numbers = input$show_numbers,
      fontsize_number = 8,
      cutree_rows = input$cutree_rows,
      cutree_cols = input$cutree_cols,
      fontsize = 18
    )
    
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })
}

# ---- Run App ----
shinyApp(ui = ui, server = server)

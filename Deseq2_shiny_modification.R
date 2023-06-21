#### Generate the count matrix data

# Generate count data
set.seed(123)
numGenes <- 1000  # Number of genes
numSamples <- 6   # Number of samples

# Create a matrix of random counts
countMatrix <- matrix(
  rpois(numGenes * numSamples, lambda = 10),
  nrow = numGenes,
  ncol = numSamples,
  dimnames = list(paste0("Gene", 1:numGenes), paste0("Sample", 1:numSamples))
)

# Create sample information
sampleInfo <- data.frame(
  condition = c(rep("Control", 3), rep("Treatment", 3)),
  replicate = rep(c("A", "B", "C"), 2)
)

# Save count data and coldata as CSV files
write.csv(countMatrix, file = "count_data.csv", row.names = TRUE)
write.csv(sampleInfo, file = "coldata.csv", row.names = FALSE)


### Shiny app

library(shiny)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Define the UI
ui <- fluidPage(
  titlePanel("DESeq2 Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("countData", "Upload Count Data (CSV format)", accept = ".csv"), # accept also txt file 
      fileInput("colData", "Upload ColData (CSV format)", accept = ".csv"),
      actionButton("runAnalysis", "Run DESeq2 Analysis"),
      downloadButton("downloadResults", "Download DESeq2 Results")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotOutput("heatmapPlot")),
        tabPanel("Enhanced Volcano", plotOutput("volcanoPlot")),
        tabPanel("MA Plot", plotOutput("maPlot"))
      )
    )
  )
)

# Define the server
server <- function(input, output) {
  
  # Read count data and coldata
  countData <- reactive({
    req(input$countData)
    read.csv(input$countData$datapath, row.names = 1)
  })
  
  colData <- reactive({
    req(input$colData)
    read.csv(input$colData$datapath, row.names = NULL)
  })
  
  # Run DESeq2 analysis
  dds <- eventReactive(input$runAnalysis, {
    DESeqDataSetFromMatrix(countData = countData(), colData = colData(), design = ~ condition)
  })
  
  ddsResult <- eventReactive(input$runAnalysis, {
    DESeq(dds())
  })
  # Display DESeq2 results as table
  output$resultsTable <- renderTable({
    req(input$runAnalysis)
    results(ddsResult())
  })
  
  # Save DESeq2 results as CSV
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("deseq2_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(as.data.frame(results(ddsResults())), file, row.names = TRUE)
    }
  )
  
  # Generate plots
  output$heatmapPlot <- renderPlot({
    req(input$runAnalysis)
    counts <- count(dds, normalized = TRUE)
    pheatmap(counts, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation")
  })
  
  output$volcanoPlot <- renderPlot({
    req(input$runAnalysis)
    res <- results(ddsResults())
    enhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', xlim = c(-10, 10), ylim = c(0, 10))
  })
  
  output$maPlot <- renderPlot({
    req(input$runAnalysis)
    res <- results(ddsResults())
    ggplot(res, aes(x = log2FoldChange, y = baseMean)) +
      geom_point(alpha = 0.6, size = 1, color = "#2c7fb8") +
      theme_bw() +
      xlim(c(-10, 10)) +
      ylim(c(0, max(res$baseMean))) +
      labs(x = "log2 Fold Change", y = "Base Mean")
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

#devtools::install_github('Linastep-bot/projectMaster')
library(projectMaster)
package_manager(c('shiny', 'DESeq2', 'shinyjs', 'ggplot2', 'pheatmap', 'EnhancedVolcano', "clusterProfiler", "tidyverse", "readr", "enrichplot", "data.table"))

options(shiny.maxRequestSize=100*1024^2)

ui <- fluidPage(
  navbarPage("NO",
  tabPanel("DESeq2 from count matrix", 
  sidebarLayout(
    sidebarPanel(
      fileInput("countData", "Upload Count Data (CSV format)", accept = ".csv"),
      fileInput("colData", "Upload ColData (CSV format)", accept = ".csv"),
      actionButton("runAnalysis", "Run DESeq2 Analysis"),
      downloadButton("downloadResults", "Download DESeq2 Results")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("DESeq2 Results", tableOutput("resultsTable")),
        tabPanel("Heatmap", plotOutput("heatmapPlot")),
        tabPanel("Enhanced Volcano", plotOutput("volcanoPlot")),
        tabPanel("MA Plot", plotOutput("maPlot"))
      )
    )
    )
  ),

    tabPanel("DESeq2 from salmon files", 
    sidebarLayout(
      sidebarPanel(
        fileInput("salmon_files", "Upload Salmon Files", multiple = TRUE, accept = ".sf"),
        fileInput("metadata", "Upload Metadata", accept = c(".xlsx", ".txt", "xls", ".csv")),
        fileInput("gff", "Upload GFF", accept = c(".gff", ".gff3")),
        fileInput("fasta", "Upload Reference Genome", accept = c(".fa", ".fasta")),
        actionButton("runAnalysis2", "Run DESeq2"),
        downloadButton("downloadResultssal", "Download DESeq2 Results")
      ),
    
    mainPanel(
      mainPanel(
        tabsetPanel(
          tabPanel("Heatmap", plotOutput("heatmapPlot1")),
          tabPanel("Enhanced Volcano", plotOutput("volcanoPlot1")),
          tabPanel("MA Plot", plotOutput("maPlot1"))
        )
      ))
      )
    ),
  
  tabPanel("Cytoscape",
             fluidRow(
               column(
                 width = 12,
                 tags$iframe(src = "https://www.ndexbio.org/iquery/", width = "100%", height = "950px")
               )
             )
           
           ),
  
  
  tabPanel("KEGG",
           sidebarLayout(
             sidebarPanel(
               selectInput("FilterBy", label = "Filter By",
                           choices = c("Overexpressed genes", "Underexpressed genes"),
                           selected = "Overexpressed genes")
             ),
             mainPanel(
               h4("KEGG"),
               splitLayout(
                 plotOutput("dotplot"),
                 plotOutput("emapplot")
               ),
               fluidRow(
                 column(3, downloadButton('downloadPlot1', "Download Dotplot")),
                 column(3, offset = 3, downloadButton('downloadPlot2', "Download Emapplot"))
               )
             )
           )
  ),
  
  
  
  tabPanel("About",
           includeMarkdown("README.md")
  )
  
  
  
  )
  
  )



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
  dds <- reactive({
    req(input$runAnalysis)
    DESeqDataSetFromMatrix(countData = countData(), colData = colData(), design = ~ as.factor(colData()[,1]))
  })
  
  ddsResults <- reactive({
    DESeq(dds())
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
    counts <- counts(estimateSizeFactors(ddsResults()), normalized=TRUE)
    pheatmap(counts, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation")
  })
  
  # Display DESeq2 results as table
  output$resultsTable <- renderTable({
    req(input$runAnalysis)
    res <- results(ddsResults())
    res
  })
  
  
  output$volcanoPlot <- renderPlot({
    req(input$runAnalysis)
    res <- results(ddsResults())
    EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', xlim = c(-10, 10), ylim = c(0, 10))
  })
  
  output$maPlot <- renderPlot({
    req(input$runAnalysis)
    res <- as.data.frame(results(ddsResults()))
    ggplot(res, aes(y = log2FoldChange, x = baseMean)) +
      geom_point(alpha = 0.6, size = 1, color = "#2c7fb8") +
      theme_bw() +
      ylim(c(-10, 10)) +
      xlim(c(0, max(res$baseMean))) +
      labs(y = "log2 Fold Change", x = "Base Mean")
  })
  
  
  # Import dataset
  deseq2_results <- read_csv("deseq2_results.csv")
  colnames(deseq2_results)[colnames(deseq2_results) == "...1"] <- "ID"
  deseq2_results_2 <- deseq2_results %>% drop_na(log2FoldChange)
  
  # Create geneUniverse
  geneUniverse <- deseq2_results_2$ID
  
  observeEvent(input$FilterBy, {
    # Filter based on overexpressed or underexpressed genes
    if (input$FilterBy == "Overexpressed genes") {
      deGenes <- deseq2_results_2[deseq2_results_2$padj < 0.05 & abs(deseq2_results_2$log2FoldChange) > log2(2), ]
      deGenes <- deGenes %>% drop_na(ID)
      deGenes <- deGenes$ID
    } else {
      deGenes <- deseq2_results_2[deseq2_results_2$padj < 0.05 & abs(deseq2_results_2$log2FoldChange) < log2(2), ]
      deGenes <- deGenes %>% drop_na(ID)
      deGenes <- deGenes$ID
    }
    
    # KEEG enrichment
    ans.kegg <- enrichKEGG(gene = deGenes,
                           organism = 'sce', 
                           universe = geneUniverse,
                           pvalueCutoff = 0.05)
    tab.kegg <- as.data.frame(ans.kegg)
    tab.kegg <- subset(tab.kegg, Count > 5)
    
    # Plot KEGG results - dotplot
    output$dotplot <- renderPlot({
      dotplot(ans.kegg, showCategory = 10) + ggtitle("KEGG")
    })
    
    # Plot KEGG results - emapplot
    output$emapplot <- renderPlot({
      x2 <- pairwise_termsim(ans.kegg) 
      emapplot(x2, showCategory = 10, cex.params = list(category_node = 1, category_label = 0.6, line = 0.5))
    })
  })
  
  # Download dotplot
  output$downloadPlot1 <- downloadHandler(
    filename = function() {
      paste("dotplot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = output$dotplot(), device = "png")
    }
  )
  
  
}










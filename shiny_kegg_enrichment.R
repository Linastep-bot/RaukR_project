#remotes::install_github("Linastep-bot/projectMaster")
library(projectMaster)
package_manager(c("shiny", "clusterProfiler", "tidyverse", "readr", "enrichplot", "data.table"))


# UI
ui <- fluidPage(
  navbarPage("KEGG enrichment analysis",
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
             )
  )
)

# Server
server <- function(input, output) {
  
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
    # Create an empty vector to store the split results
    split_results <- character(length(ans.kegg@result$Description))
    
    # For loop to split and keep the first part
    for (i in seq_along(ans.kegg@result$Description)) {
      split_parts <- unlist(strsplit(ans.kegg@result$Description[i], "-"))[[1]]  # Split the string
      split_results[i] <- split_parts[1]  # Keep the first part
    }
    
    # Convert the list to a vector if desired
    ans.kegg@result$Description <- unlist(split_results)
    
    # Plot KEGG results - dotplot
    output$dotplot <- renderPlot({
      dotplot(ans.kegg, showCategory = 10) + ggtitle("KEGG")
    })
    
    # Plot KEGG results - emapplot
    output$emapplot <- renderPlot({
      x2 <- pairwise_termsim(ans.kegg) 
      emapplot(x2, showCategory = 10, cex.params = list(category_node = 1, category_label = 0.6, line = 0.5), height = "600px")
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

# Run the Shiny app
shinyApp(ui = ui, server = server)
         
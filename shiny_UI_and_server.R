#evtools::install_github('Linastep-bot/projectMaster')
library(projectMaster)
package_manager(c('shiny', 'DESeq2', 'shinyjs'))

options(shiny.maxRequestSize=100*1024^2)

ui <- fluidPage(
  navbarPage("NO",
  tabPanel("DESeq2 from count matrix", 
  sidebarLayout(
    sidebarPanel(
      fileInput("countMatrix", "Upload Count Matrix", accept = c(".csv", ".txt")),
      fileInput("designMatrix", "Upload Design Matrix", accept = c(".csv", ".txt")),
      actionButton("runAnalysis", "Run DESeq2")
    ),
    
    mainPanel(
      textOutput("results"))
    )
  ),

    tabPanel("DESeq2 from salmon files", 
    sidebarLayout(
      sidebarPanel(
        fileInput("salmon_files", "Upload Salmon Files", multiple = TRUE, accept = ".sf"),
        fileInput("metadata", "Upload Metadata", accept = c(".xlsx", ".txt", "xls", ".csv")),
        fileInput("gff", "Upload GFF", accept = c(".gff", ".gff3")),
        fileInput("fasta", "Upload Reference Genome", accept = c(".fa", ".fasta")),
        actionButton("runAnalysis2", "Run DESeq2")
      ),
    
    mainPanel(
      textOutput("results"))
      )
    ),
  
  tabPanel("Cytoscape",
             fluidRow(
               column(
                 width = 12,
                 tags$iframe(src = "https://www.ndexbio.org/iquery/", width = "100%", height = "950px")
               )
             )
           
           )
  
  )
)


server <- function(input, output, session) {
  

  observeEvent(input$countMatrix, {
    # Read the count matrix file and extract individual names
    countMatrix <- read.csv(input$countMatrix$datapath)
    individuals <- colnames(countMatrix)
    
    # Update the choices of the checkbox group
    updateCheckboxGroupInput(session, "individuals", choices = individuals)
  })
  
  
  observeEvent(input$runAnalysis, {
    # Read the count matrix and design matrix files
    countMatrix <- read.csv(input$countMatrix$datapath)
    designMatrix <- read.csv(input$designMatrix$datapath)
    
    # Extract the selected individuals
    selectedIndividuals <- input$individuals
    
    # Subset the count matrix to include only selected individuals
    countMatrixSubset <- countMatrix[, selectedIndividuals]
    
    # Perform DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = countMatrixSubset, colData = designMatrix, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    
    # Display the results
    output$results <- renderText({
      # Format and display the differential expression results
      # Customize this based on your specific needs
      paste("Differential Expression Results:")
      head(res)
    })
  })
  
  observeEvent(input$runAnalysis2, {
    # Read the count matrix and design matrix files
    countMatrix <- read.csv(input$countMatrix$datapath)
    designMatrix <- read.csv(input$designMatrix$datapath)
    
    # Extract the selected individuals
    selectedIndividuals <- input$individuals
    
    # Subset the count matrix to include only selected individuals
    countMatrixSubset <- countMatrix[, selectedIndividuals]
    
    # Perform DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = countMatrixSubset, colData = designMatrix, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    
    # Display the results
    output$results <- renderText({
      # Format and display the differential expression results
      # Customize this based on your specific needs
      paste("Differential Expression Results:")
      head(res)
    })
  })
}


shinyApp(ui = ui, server = server)




































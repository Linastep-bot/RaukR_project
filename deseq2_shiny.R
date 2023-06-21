remotes::install_github('Linastep-bot/projectMaster')
library(projectMaster)
package_manager(c("shiny", "DESeq2"))

ui <- fluidPage(
  titlePanel("DESeq2 Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("countMatrix", "Upload Count Matrix", accept = ".csv"),
      fileInput("designMatrix", "Upload Design Matrix", accept = ".csv"),
      checkboxGroupInput("individuals", "Select Individuals", choices = NULL),
      actionButton("runAnalysis", "Run DESeq2")
    ),
    mainPanel(
      textOutput("results")
    )
  )
)


server <- function(input, output) {
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
}


shinyApp(ui = ui, server = server)



#my edits
ui <- fluidPage(
  titlePanel("DESeq2 Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("salmon_files", "Upload Salmon Files", multiple = TRUE, accept = ".sf"),
      fileInput("salmon_files", "Upload Metadata", accept = c(".xlsx", ".txt", "xls", ".csv")),
      fileInput("gff", "Upload GFF", accept = c(".gff", ".gff3")),
      fileInput("fasta", "Upload Reference Genome", accept = c(".fa", ".fasta")),
      checkboxGroupInput("individuals", "Select Individuals", choices = NULL),
      actionButton("runAnalysis", "Run DESeq2")
    ),
    mainPanel(
      textOutput("results")
    )
  )
)


server <- function(input, output) {
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
    
    #Filtering matrix
    
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

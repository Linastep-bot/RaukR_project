library(Biostrings)
library(RCircos)
library(dplyr)
library(plotly)

align_contigs <- function(assembly_a_file, assembly_b_file, bin_size) {
  # Read genome assemblies
  assembly_a <- readDNAStringSet(assembly_a_file)
  assembly_b <- readDNAStringSet(assembly_b_file)
  
  # Initialize a data frame to store the best alignments
  alignment_df <- data.frame(contig_a = character(),
                             contig_b = character(),
                             score = numeric(),
                             header_a = character(),
                             header_b = character(),
                             stringsAsFactors = FALSE)
  
  # Perform sequence alignment for each contig in assembly A to find the best match in assembly B
  for (i in seq_along(assembly_a)) {
    contig_a <- assembly_a[[i]]
    best_score <- -Inf
    best_match <- ""
    header_a <- names(assembly_a)[i]
    header_b <- ""
    
    for (j in seq_along(assembly_b)) {
      contig_b <- assembly_b[[j]]
      alignment <- pairwiseAlignment(contig_a, contig_b)
      
      if (alignment@score > best_score) {
        best_score <- alignment@score
        best_match <- contig_b
        header_b <- names(assembly_b)[j]
      }
    }
    
    alignment_df <- rbind(alignment_df, data.frame(contig_a = as.character(contig_a),
                                                   contig_b = as.character(best_match),
                                                   score = best_score,
                                                   header_a = header_a,
                                                   header_b = header_b,
                                                   stringsAsFactors = FALSE))
  }
  
  return(alignment_df)
}


# Usage example
assembly_a_file <- "path/to/fasta"
assembly_b_file <- "path/to/fasta"
bin_size <- 50

alignment_result <- align_contigs(assembly_a_file, assembly_b_file, bin_size)

# Create the interactive synteny plot
synteny_plot <- plot_ly(data = alignment_result, x = ~header_a, xaxis = "x", y = ~header_b, yaxis = "y",
                        text = ~paste("Score: ", score),
                        hoverinfo = "text", mode = "markers",
                        marker = list(symbol = "line-ew-open", size = 10, color = ~score, colorscale = "Viridis"))

# Set plot layout
layout <- list(
  xaxis = list(title = "Contig A", tickangle = 45),
  yaxis = list(title = "Contig B", tickangle = 45),
  hovermode = "closest"
)

# Combine plot and layout
synteny_plot <- synteny_plot %>% layout(layout)

# Display the plot
synteny_plot

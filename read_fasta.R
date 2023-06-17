readMultiFasta <- function(file_path) {
  # Open the file for reading
  con <- file(file_path, "r")
  
  # Initialize variables
  sequences <- list()
  curr_sequence_id <- NULL
  curr_sequence <- ""
  
  # Read the file line by line
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    line <- gsub("^\\s+|\\s+$", "", line)  # Remove leading/trailing whitespace
    
    if (nchar(line) > 0) {
      if (substr(line, 1, 1) == ">") {
        # New sequence identifier encountered
        if (!is.null(curr_sequence_id)) {
          # Store previous sequence and identifier
          sequences[[curr_sequence_id]] <- curr_sequence
        }
        # Reset current sequence variables
        curr_sequence_id <- substring(line, 2)
        curr_sequence <- ""
      } else {
        # Append sequence line to the current sequence
        curr_sequence <- paste0(curr_sequence, line)
      }
    }
  }
  
  # Store the last sequence if it exists
  if (!is.null(curr_sequence_id) && nchar(curr_sequence) > 0) {
    sequences[[curr_sequence_id]] <- curr_sequence
  }
  
  # Close the file connection
  close(con)
  
  # Return the list of sequences
  return(sequences)
}

#USAGE
fasta_file <- "path/to/your/file.fasta"
sequences <- readMultiFasta(fasta_file)

str(sequences)

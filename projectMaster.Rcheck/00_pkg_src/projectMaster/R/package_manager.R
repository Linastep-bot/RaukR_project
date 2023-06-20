#install all packages and load them 

package_manager <- function(packages) {
  # Check if packages are installed
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    # Install missing packages
    
    # Check if the package name is present in the list of available packages
    cran_packages <- available.packages()
    for (pac in missing_packages) {
      if (pac %in% cran_packages[, "Package"]) {
        install.packages(pac)
      } else {
        cat("Package", pac, "is not available in CRAN. Installing from BiocManager\n")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
          library(BiocManager)
        }
        
        BiocManager::install(pac, update = TRUE, ask = FALSE)
      }
    }
  }
  
  # Load the installed packages
  for (pkg in packages) {
    library(pkg, character.only = TRUE)
  }
}


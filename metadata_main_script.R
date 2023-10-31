# Functions to export

column_usage_check <- function(dt, step_string = "first standardization") {
  a <- lapply(dt, function(x) x != "" & !is.na(x))
  setDT(a)

  colSums(a)
  paste0(round(colSums(a) / nrow(a), 2)*100, "%")
  b <- noquote(paste0(round(colSums(a) / nrow(a), 2)*100, "%"))
  names(b) <- colnames(dt)
  print(paste("Column usage at:", step_string))
  print(b)
  return(b)
}

# Checks the R session directory is the main "riboseqorg_metadata" dir
if (basename(getwd()) != "riboseqorg_metadata") {
  cat("Warning: The working directory is not set to 'riboseqorg_metadata'. 
      temp_files will be produced in current directory.\n")
} 

# set up temp file directory within project
dir.create(file.path(dirname("."), "temp_files"))
dir.create('SraRunInfo')

################################################################################
#
# If downloads have already been carried out, change "retrieve_riboseq" to FALSE
# 
################################################################################

retrieve_riboseq <- TRUE

if ( retrieve_riboseq == TRUE) {
  cat ("Finding and fetching Ribo-seq studies.\n")
# 0.1 Find and Fetch 
  source("metadata_find&fetch.R")
# 0.2 Identify Ribo-Seq
  source("finding_riboseq.R")
  
} else {
  cat("No new Ribo-seq fetched.\n")
}


# 1. Merge columns
source("metadata_cleanup_columns.R")
# 2. Automatic Standardize column values
source("metadata_standardize_column_values.R")
# 3. Half manual standardize column values
source("metadata_standardized_columns_post_cleanup.R")


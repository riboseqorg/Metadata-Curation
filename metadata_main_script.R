#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

#' Check column usage in a data.table
#' 
#' @param dt data.table to analyze
#' @param step_string String indicating current processing step
#' @return Named vector of column usage percentages
#' @export
column_usage_check <- function(dt, step_string = "first standardization") {
  if (!is.data.table(dt)) setDT(dt)
  
  # Calculate non-empty values
  usage <- dt[, lapply(.SD, function(x) mean(!is.na(x) & x != "") * 100)]
  
  # Format percentages
  percentages <- sprintf("%.2f%%", as.numeric(usage))
  names(percentages) <- names(usage)
  
  # Log results
  log_message <- sprintf("Column usage at: %s\n%s", 
                        step_string,
                        paste(names(percentages), percentages, sep = ": ", collapse = "\n"))
  message(log_message)
  
  return(percentages)
}

#' Main workflow function
#' 
#' @param output_dir Directory for output files
#' @param temp_dir Directory for temporary files
#' @param sra_dir Directory for SRA run info
#' @param retrieve_new Logical, whether to retrieve new data
#' @param script_dir Directory containing processing scripts
#' @return Invisible NULL
#' @export
main <- function(output_dir, temp_dir, sra_dir, retrieve_new = TRUE, script_dir = NULL) {
  # Create directories if they don't exist
  dirs <- c(output_dir, temp_dir, sra_dir)
  sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  
  # Set up logging
  log_file <- file.path(output_dir, "process_log.txt")
  con <- file(log_file, "w")
  sink(con, type = "message")
  on.exit({sink(type = "message"); close(con)})
  
  # Log start time and parameters
  message(sprintf("Process started at: %s", Sys.time()))
  message(sprintf("Output directory: %s", output_dir))
  message(sprintf("Temporary directory: %s", temp_dir))
  message(sprintf("SRA directory: %s", sra_dir))
  
  tryCatch({
    if (retrieve_new) {
      message("Finding and fetching Ribo-seq studies...")
      print("finding")
      
      # Source scripts with error handling
      source_script <- function(script_name) {
        script_path <- if (!is.null(script_dir)) {
          file.path(script_dir, script_name)
        } else {
          script_name
        }
        
        if (!file.exists(script_path)) {
          stop(sprintf("Script not found: %s", script_path))
        }
        
        message(sprintf("Running script: %s", script_name))
        source(script_path, local = TRUE)
      }
      # Execute pipeline steps
      print("find and fetch")
      source_script("metadata_find&fetch.R")
      print("finding_riboseq")
      source_script("finding_riboseq.R")
    } else {
      message("Skipping new Ribo-seq fetching.")
    }
    
    # Process metadata
    print("cleanup")
    source_script("metadata_cleanup_columns.R")
    print("standardise")
    source_script("metadata_standardize_column_values.R")
    print("post_cleanup")
    source_script("metadata_standardized_columns_post_cleanup.R")
    
    message(sprintf("Process completed successfully at: %s", Sys.time()))
    
  }, error = function(e) {
    message(sprintf("Error occurred: %s", e$message))
    stop(e)
  })
}

# Command line argument parsing
option_list <- list(
  make_option(c("-o", "--output"), type="character", default="output",
              help="Output directory [default=%default]"),
  make_option(c("-t", "--temp"), type="character", default="temp_files",
              help="Temporary files directory [default=%default]"),
  make_option(c("-s", "--sra"), type="character", default="SraRunInfo",
              help="SRA run info directory [default=%default]"),
  make_option(c("-r", "--retrieve"), type="logical", default=TRUE,
              help="Whether to retrieve new data [default=%default]"),
  make_option(c("-d", "--scripts"), type="character", default=NULL,
              help="Directory containing processing scripts [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list)
opts <- parse_args(opt_parser)

# Execute main function with command line arguments
main(
  output_dir = opts$output,
  temp_dir = opts$temp,
  sra_dir = opts$sra,
  retrieve_new = opts$retrieve,
  script_dir = opts$scripts
)
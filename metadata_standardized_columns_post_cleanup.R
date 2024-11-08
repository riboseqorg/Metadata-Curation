# Load required libraries
suppressPackageStartupMessages({
  library(ORFik)
  library(data.table)
  library(readxl)
  library(BiocParallel)
  library(googledrive)
  library(massiveNGSpipe)
})

# Constants
MISSING_VALUES <- c("none", "None", "Missing", "missing", "NA", "n/a", "N/A", 
                   "not applicable", "No", "not determined")
INFO_COLS <- c("sample_title", "Info", "sample_source", "LibraryName")
RESOURCE_PATH <- "/data/resources"
TEMP_FILES_PATH <- "/data/temp_files"

#' Generate file path with current date
#' @param base_path Base path for the file
#' @param prefix Prefix for the filename
#' @param suffix Suffix for the filename
#' @return Complete file path with date
generate_dated_filepath <- function(base_path, prefix, suffix = ".csv") {
  date_str <- format(Sys.Date(), "%Y-%m-%d")
  file.path(base_path, paste0(prefix, "_", date_str, suffix))
}

#' Validate core columns against dataset
#' @param core_cols Core columns data.table
#' @param dataset Main dataset
#' @return Validated core_cols
validate_core_columns <- function(core_cols, dataset) {
  # Validate column presence
  stopifnot(
    "Core columns missing from dataset" = 
      all(core_cols[,1][[1]] %in% colnames(dataset)),
    "Mapping columns missing from dataset" = 
      all(unlist(core_cols[,2], recursive = TRUE) %in% colnames(dataset))
  )
  
  # Remove specific columns
  core_cols[!(Column %in% c("STAGE", "GENE")),]
}

#' Process library types
#' @param data Dataset to process
#' @param content_libtype Library type content mapping
#' @return Updated dataset with processed library types
process_library_types <- function(data, content_libtype) {
  for (i in seq(nrow(content_libtype))) {
    subset <- data$LIBRARYTYPE_st == "" & 
              data$LIBRARYTYPE == content_libtype$`All Names`[i]
    data[subset, LIBRARYTYPE_st := content_libtype$`Main Name`[i]]
  }
  
  # Process specific patterns
  patterns <- list(
    list(pattern = "40s|40S", col = "Info|sample_title", value = "40S"),
    list(pattern = "[Ribo]", col = "sample_title", value = "RFP", fixed = TRUE),
    list(pattern = "Ribosome Profiling", col = "INHIBITOR", value = "RFP"),
    list(pattern = "input|Total RNA", col = "FRACTION|sample_title", value = "RNA"),
    list(pattern = "ITP|ribo_mesc", col = "sample_title", value = "RFP"),
    list(pattern = "RMS", col = "sample_title", value = "RMS"),
    list(pattern = "Snap25_ip", col = "sample_title", value = "RIP"),
    list(pattern = "RP_mRNA", col = "sample_title", value = "RNA"),
    list(pattern = "RF_M", col = "sample_title", value = "RFP")
  )
  
  for (p in patterns) {
    cols <- strsplit(p$col, "\\|")[[1]]
    for (col in cols) {
      data[LIBRARYTYPE_st == "" & 
           grepl(p$pattern, get(col), fixed = p$fixed), 
           LIBRARYTYPE_st := p$value]
    }
  }
  
  # Process bioproject specific cases
  bioproject_mapping <- list(
    "PRJNA472972" = "RFP",
    "PRJNA579539" = "LSU"
  )
  
  for (proj in names(bioproject_mapping)) {
    data[BioProject == proj & LIBRARYTYPE_st == "", 
         LIBRARYTYPE_st := bioproject_mapping[[proj]]]
  }
  
  data
}

#' Process sex information
#' @param data Dataset to process
#' @return Updated dataset with processed sex information
process_sex_information <- function(data) {
  # Clean missing values
  data[Sex %in% MISSING_VALUES, Sex := ""]
  
  # Standardize existing values
  data[Sex == "Male", Sex := "male"]
  data[Sex %in% c("pooled male and female", "mixed"), 
       Sex := "Mix of male/female"]
  
  # Set sex based on cell lines
  female_cell_lines <- c("HeLa")
  male_cell_lines <- c("Hek293")
  
  data[CELL_LINE_st %in% female_cell_lines, Sex := "female"]
  data[CELL_LINE_st %in% male_cell_lines, Sex := "male"]
  
  data
}

#' Main processing function
#' @param input_path Path to input file
#' @return Processed dataset
main <- function() {
  # Read input data
  input_path <- file.path(TEMP_FILES_PATH, "standardized_columns_with_original.csv")
  output_path <- generate_dated_filepath(TEMP_FILES_PATH, "standardized_columns_final")
  
  # Read main dataset
  x_st <- fread(input_path)
  # Read and process core columns
  core_cols <- fread(file.path(RESOURCE_PATH, "Column_value_mapping.csv"),
                    colClasses = c("character", "character"))

  setDT(core_cols)
  core_cols[, 2] <- lapply(seq(nrow(core_cols)), function(y) 
    strsplit(core_cols[y,2][[1]][[1]], split = ", ")[[1]])
  
  core_cols <- validate_core_columns(core_cols, x_st)
  
  # Read content mapping
  content <- fread(file.path(RESOURCE_PATH, "Content.csv"),
                  colClasses = rep("character", 4))
  
  # Process library types
  content_libtype <- content[content$Column == "LIBRARYTYPE",]
  x_final <- copy(x_st)

  # Add SRA IDs
  ids <- fread(file.path(TEMP_FILES_PATH, "SRA_ids.csv"))
  x_final <- cbind(x_final, ids)
  
  # Process data
  x_final <- process_library_types(x_final, content_libtype)
  x_final <- process_sex_information(x_final)
  
  # Generate unique names and check duplicates
  x_final[, name := do.call(paste, c(.SD, sep = "_")), 
          .SDcols = grep("_st", colnames(x_final))]
  x_final[, not_unique := duplicated(name), 
          by = .(BioProject, ScientificName)]
  
  # Save results
  fwrite(x_final, output_path)
  # Return statistics
  list(
    missing_library_types = nrow(x_final[LIBRARYTYPE_st == "",]),
    total_rows = nrow(x_final),
    duplicate_count = sum(x_final$not_unique)
  )
}

# Execute main function if running as script
main()

# Load required libraries
library(ORFik)
library(data.table)
library(readxl)

# Read initial metadata
metadata_df <- fread("/data/temp_files/RiboSeq_Metadata_All_Columns.csv")

# Define standard value mappings
seq_type_mappings <- list(
  "ribosome-protected" = c(Yes = "Ribo-seq", No = "mRNA-seq"),
  "depolarised" = c(Yes = "depolarised")
)

# Convert sequencing type indicators
metadata_df[`ribosome-protected` == "Yes", `ribosome-protected` := "Ribo-seq"]
metadata_df[`ribosome-protected` == "No", `ribosome-protected` := "mRNA-seq"]
metadata_df[`depolarised` == "Yes", `depolarised` := "depolarised"]

# Standardize time point values
valid_timepoints <- c(
  "1", "2", "4", "6", "8", "5", "2.5"
)
for (timepoint in valid_timepoints) {
  column_name <- paste0("Experimental Factor: ", timepoint)
  metadata_df[get(column_name) == "time point", (column_name) := timepoint]
}

# Standardize experimental factors
metadata_df[`Experimental Factor: MHV-A59` == "infect", `Experimental Factor: MHV-A59` := "MHV-A59"]
metadata_df[`scriptminer index` != "", `scriptminer index` := paste("scriptminer: index", `scriptminer index`)]

# Process virus experimental factors
virus_columns <- grep("Experimental Factor: encephalomyocarditis virus", colnames(metadata_df), value = F)
for (col_index in virus_columns) {
  current_col <- metadata_df[,col_index, with = F][[1]]
  metadata_df[current_col != "", which(seq_along(colnames(metadata_df)) == col_index) := gsub(".*\\(|\\)", "", colnames(metadata_df)[col_index])]
}

# Remove empty columns
remove_empty_columns <- function(df) {
  non_empty_counts <- sapply(df, function(col) sum(!(is.na(col) | as.character(col) == "")))
  df[, which(non_empty_counts > 0), with = F]
}
metadata_df <- remove_empty_columns(metadata_df)

# Handle duplicate columns
duplicate_columns <- table(colnames(metadata_df)) > 1
duplicate_column_names <- names(duplicate_columns[duplicate_columns])

filtered_data <- copy(metadata_df)

for (i in duplicate_column_names) {
  hits <- which(colnames(filtered_data) == i)
  filtered_data[, key_ := do.call(paste, c(.SD, sep = "_")), .SDcols = hits]
  filtered_data[, key_ := gsub("_$", "", key_)]
  filtered_data[, which(colnames(filtered_data) == i):=NULL]
  colnames(filtered_data)[which(colnames(filtered_data) == "key_")] <- i
}

stopifnot(sum(table(colnames(filtered_data)) > 1) == 0)

# Define sheet paths for rules processing
sheet_files <- list(
  core = "/data/resources/Core.csv",
  open = "/data/resources/Open.csv",
  technical = "/data/resources/Technical.csv"
)

# Process each sheet's rules
for (s in names(sheet_files)) {
  col_collapse_rules <- fread(sheet_files[[s]])
  # Process rules based on column count
if (ncol(col_collapse_rules) == 4) {
    # Core sheet (has "category" to remove, also 2 search columns)
    dt <- as.data.table(col_collapse_rules)[, 2:4]
    dt[,2] <- lapply(seq(nrow(dt)), function(y) {
      values <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T),
             gsub(pattern = '"', "", substr(strsplit(dt[y,3][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T))
      values <- values[!is.na(values)]
    })
  } else {
    # Other sheets
    dt <- as.data.table(col_collapse_rules)
    dt[,2] <- lapply(seq(nrow(dt)), function(y) {
      values <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T))
      values <- values[!is.na(values)]
    })
  }

  # Apply rules to merge columns
  for (i in seq(nrow(dt))) {
    ref <- unlist(dt[i,2], recursive = T, use.names = F)
    hits <- which(colnames(filtered_data) %in% ref)
    filtered_data[, key_ := do.call(paste, c(.SD, sep ="_")), .SDcols = hits]
    
    # Clean up merged values
    filtered_data[, key_ := gsub("^NA_|_NA$|^NA$|_NA_|_$|NA_NA$|__", "", key_)]
    filtered_data[, key_ := gsub("^_|_$", "", key_)]
    filtered_data[, key_ := gsub("^NANA$|^NA_$|^NA$|^_$", "", key_)]
    filtered_data[, key_ := gsub("_NA$|_$", "", key_)]
    
    # Remove original columns and rename merged column
    filtered_data[, which(colnames(filtered_data) %in% ref):=NULL]
    colnames(filtered_data)[which(colnames(filtered_data) == "key_")] <- dt[i,1][[1]]
  }
}

# Clean REPLICATE and TIMEPOINT columns
for (col in c("REPLICATE", "TIMEPOINT")) {
  filtered_data[, (col) := gsub("^_|_$", "", get(col))]
  filtered_data[, (col) := gsub("^NA|NA$", "", get(col))]
  filtered_data[, (col) := gsub("^_|_$", "", get(col))]
}

# Remove irrelevant columns
irrelevant <- fread("/data/resources/Irrelevant.csv")$Delete
filtered_data[, which(colnames(filtered_data) %in% irrelevant):=NULL]

# Remove duplicated Sample identifiers
identifier_patterns <- "^ENA |^ENA-|^INSDC |date|GEO Accession|Experiment Date"
filtered_data <- filtered_data[, -grep(identifier_patterns, colnames(filtered_data)), with = F]

# Save processed data
fwrite(filtered_data, file = "/data/temp_files/filtered_riboseq_done.csv")

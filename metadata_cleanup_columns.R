# Load required libraries
library(ORFik)
library(data.table)
library(readxl)

# Initialize tracking of column counts
total_columns <- c()

# Read initial data
a <- fread("/data/temp_files/RiboSeq_Metadata_All_Columns.csv")
dim(a)
total_columns <- c(total_columns, dim(a)[2])

# Convert values to meaningful strings
boolean_conversions <- list(
  "ribosome-protected" = c(Yes = "Ribo-seq", No = "mRNA-seq"),
  "depolarised" = c(Yes = "depolarised")
)
a[`ribosome-protected` == "Yes", `ribosome-protected` := "Ribo-seq"]
a[`ribosome-protected` == "No", `ribosome-protected` := "mRNA-seq"]
a[`depolarised` == "Yes", `depolarised` := "depolarised"]

# Convert time points
timepoint_conversions <- c(
  "1", "2", "4", "6", "8", "5", "2.5"
)
for (tp in timepoint_conversions) {
  col <- paste0("Experimental Factor: ", tp)
  a[get(col) == "time point", (col) := tp]
}

# Handle special cases
a[`Experimental Factor: MHV-A59` == "infect", `Experimental Factor: MHV-A59` := "MHV-A59"]
a[`scriptminer index` != "", `scriptminer index` := paste("scriptminer: index", `scriptminer index`)]

# Process encephalomyocarditis virus columns
cols <- grep("Experimental Factor: encephalomyocarditis virus", colnames(a), value = F)
for (x in cols) {
  col <- a[,x, with = F][[1]]
  a[col != "", which(seq_along(colnames(a)) == x) := gsub(".*\\(|\\)", "", colnames(a)[x])]
}

# Remove empty columns
remove.empty.cols <- function(a) {
  columns_breakdown <- sapply(a, function(x) sum(!(is.na(x) | as.character(x) == "")))
  a[, which(columns_breakdown > 0), with = F]
}
a <- remove.empty.cols(a)
total_columns <- c(total_columns, dim(a)[2])
# Merge identical named columns
duplicated_col_names <- table(colnames(a)) > 1

duplicated_col_names <- names(duplicated_col_names[duplicated_col_names])
x <- copy(a)

for (i in duplicated_col_names) {
  hits <- which(colnames(x) == i)
  x[, key_ := do.call(paste, c(.SD, sep = "_")), .SDcols = hits]
  x[, key_ := gsub("_$", "", key_)]
  x[, which(colnames(x) == i):=NULL]
  colnames(x)[which(colnames(x) == "key_")] <- i
}
dim(x)
total_columns <- c(total_columns, dim(x)[2])
stopifnot(sum(table(colnames(x)) > 1) == 0)

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
      a <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T),
             gsub(pattern = '"', "", substr(strsplit(dt[y,3][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T))
      a <- a[!is.na(a)]
    })
  } else {
    # Other sheets
    dt <- as.data.table(col_collapse_rules)
    dt[,2] <- lapply(seq(nrow(dt)), function(y) {
      a <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T))
      a <- a[!is.na(a)]
    })
  }

  # Apply rules to merge columns
  for (i in seq(nrow(dt))) {
    ref <- unlist(dt[i,2], recursive = T, use.names = F)
    hits <- which(colnames(x) %in% ref)
    x[, key_ := do.call(paste, c(.SD, sep ="_")), .SDcols = hits]
    
    # Clean up merged values
    x[, key_ := gsub("^NA_|_NA$|^NA$|_NA_|_$|NA_NA$|__", "", key_)]
    x[, key_ := gsub("^_|_$", "", key_)]
    x[, key_ := gsub("^NANA$|^NA_$|^NA$|^_$", "", key_)]
    x[, key_ := gsub("_NA$|_$", "", key_)]
    
    # Remove original columns and rename merged column
    x[, which(colnames(x) %in% ref):=NULL]
    colnames(x)[which(colnames(x) == "key_")] <- dt[i,1][[1]]
  }
  dim(x)
  total_columns <- c(total_columns, dim(x)[2])
}

# Clean REPLICATE and TIMEPOINT columns
for (col in c("REPLICATE", "TIMEPOINT")) {
  x[, (col) := gsub("^_|_$", "", get(col))]
  x[, (col) := gsub("^NA|NA$", "", get(col))]
  x[, (col) := gsub("^_|_$", "", get(col))]
}

# Remove irrelevant columns
irrelevant <- fread("/data/resources/Irrelevant.csv")$Delete
x[, which(colnames(x) %in% irrelevant):=NULL]


total_columns <- c(total_columns, dim(x)[2])


# Remove duplicated Sample identifiers
identifier_patterns <- "^ENA |^ENA-|^INSDC |date|GEO Accession|Experiment Date"
x <- x[, -grep(identifier_patterns, colnames(x)), with = F]

total_columns <- c(total_columns, dim(x)[2])

# Name the column counts for each stage
names(total_columns) <- c("raw", "empty_filtered", "equal merged", 
                         "similar merged", "open merged", "technical merged", 
                         "irrelevant filtered", "duplicated accession")

total_columns

# Save processed data
fwrite(x, file = "/data/temp_files/filtered_riboseq_done_260623.csv")
# Optional: Create subset without SRA columns for inspection
SRA_cols <- fread("/data/resources/SRA.csv")
y <- x[,!(colnames(x) %in% names(SRA_cols)), with = F]
dim(y)

# Sort columns for inspection
sort(colnames(y))
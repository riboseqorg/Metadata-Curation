library(ORFik); library(data.table); library(readxl)

total_columns <- c()

# Update this path, all else should just work
a <- fread("temp_files/RiboSeq_Metadata_All_Columns.csv")
dim(a); total_columns <- c(total_columns, dim(a)[2])


# Convert booleans to strings
a[`ribosome-protected` == "Yes", `ribosome-protected` := "Ribo-seq"]
a[`ribosome-protected` == "No", `ribosome-protected` := "mRNA-seq"]
a[`Experimental Factor: 1` == "time point", `Experimental Factor: 1` := "1"]
a[`Experimental Factor: 2` == "time point", `Experimental Factor: 2` := "2"]
a[`Experimental Factor: 4` == "time point", `Experimental Factor: 4` := "4"]
a[`Experimental Factor: 6` == "time point", `Experimental Factor: 6` := "6"]
a[`Experimental Factor: 8` == "time point", `Experimental Factor: 1` := "8"]
a[`Experimental Factor: 5` == "time point", `Experimental Factor: 5` := "5"]
a[`Experimental Factor: 2.5` == "time point", `Experimental Factor: 2.5` := "2.5"]
a[`Experimental Factor: MHV-A59` == "infect", `Experimental Factor: MHV-A59` := "MHV-A59"]
a[`depolarised` == "Yes", `depolarised` := "depolarised"]

a[`scriptminer index` != "", `scriptminer index` := paste("scriptminer: index", `scriptminer index`)]

cols <- grep("Experimental Factor: encephalomyocarditis virus", colnames(a), value = F)

for (x in cols) {
  col <- a[,x, with = F][[1]]
  a[col != "", which(seq_along(colnames(a)) == x) := gsub(".*\\(|\\)", "", colnames(a)[x])]
}


# Remove all empty columns
remove.empty.cols <- function(a) {
  columns_breakdown <- sapply(a, function(x) sum(!(is.na(x) | as.character(x) == "")))
  a[, which(columns_breakdown > 0), with = F]
}

a <- remove.empty.cols(a)

dim(a); total_columns <- c(total_columns, dim(a)[2])


# Remove total identical named columns
duplicated_col_names <- table(colnames(a)) > 1
sum(table(colnames(a)) > 1)
duplicated_col_names <- names(duplicated_col_names[duplicated_col_names])
x <- copy(a)

for (i in duplicated_col_names) {
  hits <- which(colnames(x) == i)
  x[, key_ := do.call(paste, c(.SD, sep = "_")), .SDcols = hits]
  x[, key_ := gsub("_$", "", key_)]
  x[, which(colnames(x) == i):=NULL]
  colnames(x)[which(colnames(x) == "key_")] <- i
}

dim(x); total_columns <- c(total_columns, dim(x)[2])

stopifnot(sum(table(colnames(x)) > 1) == 0)


# Remove space/_ variant naming
sheets <- c(2,3,4) # Fixed columns

for (s in sheets) {
  print(paste("Sheet:", s))
  col_collapse_rules <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                                  sheet = s)
  col_collapse_rules <- data.table(col_collapse_rules)
  if (ncol(col_collapse_rules) == 4) {
    # Sheet 2 (has "category" to remove, also 2 search columns)
    dt <- as.data.table(col_collapse_rules)[, 2:4]
    dt[,2] <- lapply(seq(nrow(dt)), function(y) {a <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T),
                                                        gsub(pattern = '"', "", substr(strsplit(dt[y,3][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T)); a <-a[!is.na(a)]})
  } else {
    # For other sheets
    dt <- as.data.table(col_collapse_rules)
    dt[,2] <- lapply(seq(nrow(dt)), function(y) {a <- c(gsub(pattern = '"', "", substr(strsplit(dt[y,2][[1]], split = ".\\, ")[[1]], 2, 1000), fixed = T)); a <-a[!is.na(a)]})
  }

  for (i in seq(nrow(dt))) {
    ref <- unlist(dt[i,2], recursive = T, use.names = F)
    hits <- which(colnames(x) %in% ref)
    x[, key_ := do.call(paste, c(.SD, sep ="_")), .SDcols = hits]
    x[, key_ := gsub("^NA_", "", key_)]
    x[, key_ := gsub("^NA_", "", key_)]
    x[, key_ := gsub("^NA_", "", key_)]
    x[, key_ := gsub("_NA$", "", key_)]
    x[, key_ := gsub("_NA$", "", key_)]
    x[, key_ := gsub("_NA$", "", key_)]
    x[, key_ := gsub("_NA$", "", key_)]
    x[, key_ := gsub("^NA$", "_", key_)]
    x[, key_ := gsub("_NA_", "_", key_)]
    x[, key_ := gsub("_$", "", key_)]
    x[, key_ := gsub("NA_NA$", "", key_)]
    x[, key_ := gsub("_$", "", key_)]
    x[, key_ := gsub("__", "", key_)]
    x[, key_ := gsub("^_", "", key_)]
    x[, key_ := gsub("^_", "", key_)]
    x[, key_ := gsub("^NANA$", "_", key_)]
    x[, key_ := gsub("^NA_$", "_", key_)]
    x[, key_ := gsub("^NA$", "_", key_)]
    x[, key_ := gsub("^_$", "", key_)]
    x[, key_ := gsub("_$", "", key_)]
    x[, key_ := gsub("_$", "", key_)]
    x[, key_ := gsub("_NA$", "_", key_)]
    x[, which(colnames(x) %in% ref):=NULL]
    colnames(x)[which(colnames(x) == "key_")] <- dt[i,1][[1]]
  }
  dim(x); total_columns <- c(total_columns, dim(x)[2])
}

x[, REPLICATE := gsub("^_|_$", "", REPLICATE)]
x[, REPLICATE := gsub("^NA|NA$", "", REPLICATE)]
x[, REPLICATE := gsub("^_|_$", "", REPLICATE)]
x[, TIMEPOINT := gsub("^_|_$", "", TIMEPOINT)]
x[, TIMEPOINT := gsub("^NA|NA$", "", TIMEPOINT)]
x[, TIMEPOINT := gsub("^_|_$", "", TIMEPOINT)]


# Remove irrelevant / duplicates
irrelevant <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                                sheet = 5)$Delete
x[, which(colnames(x) %in% irrelevant):=NULL]
dim(x); total_columns <- c(total_columns, dim(x)[2])


# Remove duplicated Sample identifiers
x <- x[, -grep("^ENA |^ENA-|^INSDC |date|GEO Accession|Experiment Date", colnames(x)), with = F]
dim(x); total_columns <- c(total_columns, dim(x)[2])

names(total_columns) <- c("raw", "empty_filtered", "equal merged", "similar merged", "open merged", "technical merged", "irrelevant filtered", "duplicated accession")
total_columns

fwrite(x, file = "temp_files/filtered_riboseq_done_260623.csv")


# For inspection, ignore SRA columns and made columns
SRA_cols <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                                sheet = 6)
y <- x[,!(colnames(x) %in% names(SRA_cols)), with = F]
dim(y)


# # Analyse remaining with these lines:
sort(colnames(y))


# unique(z$`crosslinking`)
# View(remove.empty.cols(z[`Library Pool ID` != "",]))

# table(z$`crosslinking`)
# grep("dis", colnames(x), value = T)
# View(remove.empty.cols(x[`BioProject` == "PRJEB17636",]))
# View(remove.empty.cols(x[`crosslinking` != "",]))
# View(remove.empty.cols(z[`Library Pool ID` != "",]))
# print(dim(y))


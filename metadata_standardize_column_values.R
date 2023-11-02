library(ORFik); library(data.table); library(readxl); library(BiocParallel)
library(massiveNGSpipe)

# devtools::load_all()
findFromPath <- ORFik:::findFromPath
content <-  googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                      sheet = 1, col_types = "cccL")
content <- as.data.table(content)
colnames(content) <- c("Category","Column", "mainName", "allNames")

x <- fread("temp_files/filtered_riboseq_done_260623.csv")

# Table info
# Total rows
dim(x)
# Top 10 organisms
head(sort(table(x$ScientificName), decreasing = T), 10)
# Sequencing platforms
table(x$Platform)
# Sequencing strategies
table(x$LibraryStrategy)
# Filter non accepted rows
x <- massiveNGSpipe:::pipeline_metadata_filter(x)

# Remove SRA cols not needed now
SRA_cols <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                      sheet = 6)
SRA_cols <- colnames(SRA_cols)
SRA_cols <- SRA_cols[!(SRA_cols %in% c("sample_title", "Info", "sample_source", "LibraryName", "ScientificName"))]
SRA_to_id <- x[, c("Run", "BioProject"), with = F]
fwrite(SRA_to_id, file = "temp_files/SRA_ids.csv")
x <- x[,!(colnames(x) %in% SRA_cols), with = F]


############# Standardize columns
info_cols <- c("sample_title", "Info", "sample_source", "LibraryName")
standardize_col <- function(x, cols_to_check, check_table) {
  stopifnot(all(cols_to_check %in% colnames(x)))

  hits <- bplapply(cols_to_check, function(y) {
    print(y)
    findFromPath(x[, y, with = F][[1]], check_table, "auto")
  })
  names(hits) <- cols_to_check
  setDT(hits)
  return(hits)
}

core_cols <-  googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                      sheet = 7, col_types = "cL")
setDT(core_cols)
core_cols[,2] <- lapply(seq(nrow(core_cols)), function(y) strsplit(core_cols[y,2][[1]][[1]], split = ", ")[[1]])
stopifnot(all(core_cols[,1][[1]] %in% colnames(x)))
stopifnot(all(unlist(core_cols[,2], recursive = T) %in% colnames(x)))
core_cols <- core_cols[!(Column %in% c("STAGE", "GENE")),]

mapping_functions <- list(CELL_LINE = ORFik:::cellLineNames(),
                          TISSUE = list(ORFik:::tissueNames(), ORFik:::cellLineNames(T)),
                          INHIBITOR = ORFik:::inhibitorNames(),
                          TIMEPOINT = ORFik:::stageNames(),
                          FRACTION = ORFik:::fractionNames(),
                          REPLICATE = ORFik:::repNames(),
                          CONDITION = ORFik:::conditionNames(),
                          LIBRARYTYPE = ORFik:::libNames())
is_nested <- !unlist(lapply(mapping_functions, is.data.table))
l <- lapply(seq(nrow(core_cols)), function(y) {
  print(core_cols[y,]$Column)
  if (is_nested[y]) {
    maps <- mapping_functions[[y]]
  } else maps <- mapping_functions[y]
  lapply(maps, function(map_fun) standardize_col(x, c(core_cols[y,]$Mapping[[1]], info_cols), map_fun))
})

l[[2]] <- list(TISSUE = as.data.table(dplyr::bind_cols(l[2])))

dt <- data.table()
for (j in seq_along(l)) {
  for (i in seq_along(l[j][[1]][[1]])) {
    if (i == 1) {
      col_res <- l[j][[1]][[1]][,1][[1]]
    } else col_res[col_res == ""] <- l[j][[1]][[1]][,i, with = F][[1]][col_res == ""]
    if (all(col_res != "")) break
  }
  dt <- data.table(dt, col_res)
}
colnames(dt) <- core_cols$Column
lapply(dt, function(x) table(x))

# This function is not present in the script!
# column_usage_check(dt)

dt_st <- copy(dt)
colnames(dt_st) <- paste0(colnames(dt_st), "_st")
x_st <- cbind(x, dt_st)
fwrite(x_st, "temp_files/standardized_columns_with_original.csv")

#
# # Then check strategy
# filtered_RFP[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LIBRARYTYPE <-
#   ORFik:::findFromPath(filtered_RFP[LIBRARYTYPE == "" & !(LibraryStrategy %in%  c("RNA-Seq", "OTHER")),]$LibraryStrategy,
#                        ORFik:::libNames(), "auto")
#
# filtered_RFP$BATCH <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                            ORFik:::batchNames(), "auto")
# filtered_RFP$REPLICATE <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                                ORFik:::repNames(), "auto")
# filtered_RFP$TIMEPOINT <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                                ORFik:::stageNames(), "auto")
# filtered_RFP$TISSUE <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                             ORFik:::tissueNames(), "auto")
# # browser()
# filtered_RFP[TISSUE == "",]$TISSUE <- ORFik:::findFromPath(filtered_RFP[TISSUE == "",]$sample_source,
#                                                            ORFik:::tissueNames(), "auto")
# filtered_RFP$CELL_LINE <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                                ORFik:::cellLineNames(), "auto")
# filtered_RFP[CELL_LINE == "",]$CELL_LINE <- ORFik:::findFromPath(filtered_RFP[CELL_LINE == "",]$sample_source,
#                                                                  ORFik:::cellLineNames(), "auto")
# filtered_RFP[(TISSUE == "") & CELL_LINE != "",]$TISSUE <- ORFik:::findFromPath(filtered_RFP[(TISSUE == "") & CELL_LINE != "",]$CELL_LINE,
#                                                                                ORFik:::cellLineNames(TRUE), "auto")
# # Correct replicate names
# duplicated_samples <- any(duplicated(filtered_RFP$sample_title) & !is.na(filtered_RFP$sample_title))
# if (duplicated_samples) {
#   filtered_RFP[, has_dups := any(duplicated(sample_title) | all(REPLICATE == "")), by= .(study_accession)]
#   filtered_RFP[has_dups == TRUE, REPLICATE := seq(.N), by= .(study_accession, sample_title)]
#   filtered_RFP$has_dups <- NULL
# }
# # Knock outs
# filtered_RFP$GENE <- ""
# combs <- expand.grid(c("Δ", "KO"), c(" .*", "_.*", "$"))
# KO_combs <- paste(paste0(combs[,1], combs[,2]), collapse = "|")
# KO_hits <- grep(KO_combs, filtered_RFP$sample_title)
# filtered_RFP[KO_hits,]$CONDITION <- "KO"
# filtered_RFP[KO_hits,]$GENE <- gsub(KO_combs, "", filtered_RFP[KO_hits,]$sample_title)
# # Exclude KO hits column
# filtered_RFP$CONDITION[-KO_hits] <- ORFik:::findFromPath(filtered_RFP$sample_title[-KO_hits],
#                                                          ORFik:::conditionNames(), "auto")
# filtered_RFP$FRACTION <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                               ORFik:::fractionNames(), "auto")
# filtered_RFP$INHIBITOR <- ORFik:::findFromPath(filtered_RFP$sample_title,
#                                                ORFik:::inhibitorNames(), "auto")
# filtered_RFP[LIBRARYTYPE == "" & INHIBITOR != "", LIBRARYTYPE := "RFP"]
# # Add columns needed for checks
# filtered_RFP[, `:=`(KEEP = "", UNIQUE = "", CHECKED = "", name = "", not_unique = "")]
# cat("Studies kept:", "\n")
# cat(length(unique(filtered_RFP$Submission)), "\n")
# cat("Library types detected", "\n")
# print(table(filtered_RFP$LIBRARYTYPE)) # It can't find all (bad info)
# filtered_RFP[]
#
#
# ###################### Manual hard set
# # list(RFP = c("Ribosome Protected mRNA", "Ribosome protected mRNA", "ribosome-bound mRNA", "RNA (ribosome protected)"),
# #      SSU = c("40S"),
# #      LSU = c("80S"))

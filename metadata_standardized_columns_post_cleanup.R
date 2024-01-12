library(ORFik); library(data.table); library(readxl); library(BiocParallel)
library(googledrive); library(massiveNGSpipe)
# devtools::load_all()
findFromPath <- ORFik:::findFromPath
input_file_path <- "temp_files/standardized_columns_with_original.csv"

today_date <- format(Sys.Date(), "%Y-%m-%d")
final_file_path <- paste0("temp_files/standardized_columns_final_", today_date, ".csv")


x_st <- fread(input_file_path)

core_cols <-  googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                        sheet = 7, col_types = "cL")
setDT(core_cols)
core_cols[,2] <- lapply(seq(nrow(core_cols)), function(y) strsplit(core_cols[y,2][[1]][[1]], split = ", ")[[1]])
stopifnot(all(core_cols[,1][[1]] %in% colnames(x_st)))
stopifnot(all(unlist(core_cols[,2], recursive = T) %in% colnames(x_st)))
core_cols <- core_cols[!(Column %in% c("STAGE", "GENE")),]
info_cols <- c("sample_title", "Info", "sample_source", "LibraryName")

# Analyse miss
for (i in seq_along(core_cols$Column)) {
  col <- core_cols$Column[i]
  print(col)
  subset <- x_st[, paste0(col, "_st"), with = F][[1]] == "" & x_st[, col, with = F][[1]] != ""
  missing <- x_st[subset, c(info_cols, core_cols[Column == col,]$Mapping[[1]], paste0(col, "_st")), with = F][, col, with = F][[1]]
  print(unique(missing))
}

content <-  googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1TVjXdpyAMJBex-OWyfEdZf_nfF79bPd_T6TpXS1LnFM/edit#gid=515852084",
                                      sheet = 9, col_types = "cccc")
# Library type
# Hard coded
content_libtype <- content[content$Column == "LIBRARYTYPE",]
for (i in seq(nrow(content_libtype))) {
  subset <- x_st$LIBRARYTYPE_st == "" & x_st$LIBRARYTYPE == content_libtype$`All Names`[i]
  x_st[subset, LIBRARYTYPE_st := content_libtype$`Main Name`[i]]
}
# x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("ITP", sample_title),]
# Inhibitor
content_inhib <- content[content$Column == "INHIBITOR",]
# for (i in seq(nrow(content_libtype))) {
#   subset <- x_st$INHIBITOR_st == "" & x_st$INHIBITOR == content_libtype$`All Names`[i]
#   x_st[subset, INHIBITOR_st := content_libtype$`Main Name`[i]]
# }




# x_st[ , new := do.call(paste, c(.SD, sep = "_")), .SDcols=-which(colnames(x_st) %in% c("REPLICATE", "REPLICATE_st"))]
ids <- fread("temp_files/SRA_ids.csv")
x_final <- copy(x_st)
x_final <- cbind(x_final, ids)
# Grep coded
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("40s", Info), LIBRARYTYPE_st := "40S"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("40S", sample_title), LIBRARYTYPE_st := "40S"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("[Ribo]", sample_title, fixed = T), LIBRARYTYPE_st := "RFP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("Ribosome Profiling", INHIBITOR), LIBRARYTYPE_st := "RFP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("input|Total RNA", FRACTION), LIBRARYTYPE_st := "RNA"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("input|Total RNA", sample_title),]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("ITP", sample_title), LIBRARYTYPE_st := "RFP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("ribo_mesc", sample_title), LIBRARYTYPE_st := "RFP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("RMS", sample_title), LIBRARYTYPE_st := "RMS"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("Snap25_ip", sample_title), LIBRARYTYPE_st := "RIP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("RP_mRNA", sample_title), LIBRARYTYPE_st := "RNA"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("RF_M", sample_title), LIBRARYTYPE_st := "RFP"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("_RF_", sample_title) & seq(nrow(x_final)) %in% grep("_input_", sample_title), LIBRARYTYPE_st := "RNA"]

# Bioproject specific
x_final[BioProject == "PRJNA472972", LIBRARYTYPE_st := "RFP"]
x_final[BioProject == "PRJNA579539" & LIBRARYTYPE_st == "", LIBRARYTYPE_st := "LSU"]
x_final[LIBRARYTYPE_st == "" & BioProject == "PRJNA501872" & seq(nrow(x_final)) %in% grep("input", sample_title), LIBRARYTYPE_st := "RNA"]
x_final[LIBRARYTYPE_st == "" & seq(nrow(x_final)) %in% grep("-input-", sample_title), LIBRARYTYPE_st := "RNA"]

# Inhibitor
x_final[INHIBITOR_st == "" & seq(nrow(x_final)) %in% grep("Har|har|HRT", INHIBITOR, ) & seq(nrow(x_final)) %in% grep("cyclo|Cyclo|CHX|CYH", INHIBITOR), INHIBITOR_st := "chx_harr"]
table(x_final$INHIBITOR_st)
# Sex (TODO: Infer from cell lines)
missing_empty <- c("none", "None", "Missing", "missing", "NA", "n/a", "N/A", "not applicable","No", "not determined")
x_final[Sex %in% missing_empty, Sex := ""]
x_final[Sex %in% "Male", Sex := "male"]
x_final[Sex %in% c("pooled male and female", "mixed"), Sex := "Mix of male/female"]
# Female cell lines
x_final[CELL_LINE_st %in% c("HeLa"), Sex := "female"]
# Male cell lines
x_final[CELL_LINE_st %in% c("Hek293"), Sex := "male"]
table(x_final$Sex)


## Check and save
x_final[ , name := do.call(paste, c(.SD, sep = "_")), .SDcols=grep("_st", colnames(x_final))]
x_final[, not_unique := duplicated(name), by = .(BioProject, ScientificName)]
table(x_final$not_unique)
biogrouping <- paste("")
# Save results
fwrite(x_final, final_file_path)
googledrive::drive_upload(final_file_path, as_id("https://drive.google.com/drive/folders/1QLq31NOM1NgC4otN5664_Sdg6QvUlOwg"), overwrite = TRUE)

# Statistics and column usage of results
# View(x_final[not_unique == F,])
nrow(x_final[LIBRARYTYPE_st == "",]); View(x_final[LIBRARYTYPE_st == "",])
dtt <- x_final[, grep("_st", colnames(x_final)), with = F]
colnames(dtt) <- gsub("_st", "", colnames(dtt))
column_usage_check(dtt, "Final standardization")

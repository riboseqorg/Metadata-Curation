library(ORFik)
library(data.table)
library(dplyr)



bioprojects <- ORFik::get_bioproject_candidates("((Ribosomal footprinting) OR (Ribosome footprinting) OR (Ribosome profiling) OR ribo-seq)",add_study_title = TRUE)

# Path to whitelist CSV - you can change this path as needed
whitelisted_path <- "/usr/local/share/riboseq/resources/whitelisted_bioprojects.csv"

whitelisted_bioprojects <- fread(whitelisted_path)


# Identify rows in whitelisted_bioprojects that are not in bioprojects and add them to bioprojects
rows_to_add <- anti_join(whitelisted_bioprojects, bioprojects, by = "id")
bioprojects <- bind_rows(bioprojects, rows_to_add)

processed_records <- system("ls /data/SraRunInfo/", intern = TRUE) %>% sub("SraRunInfo_","",.) %>% sub(".csv","",.)
unprocessed <- bioprojects$id[!(bioprojects$id %in% processed_records)]
bioprojects$abstract <- ""

writeSRA <- function(candidate) {
  sra <- try(download.SRA.metadata(candidate, outdir = "/data/SraRunInfo", rich.format = TRUE))
  return(candidate)
}

rna_seq_meta <- lapply(unprocessed, function(x) writeSRA(x))

abstracts <- system("ls /data/SraRunInfo/abstract_*", intern = TRUE)
abstracts <- sapply(abstracts, function(x) gsub("\"","",readLines(x)[[2]]), USE.NAMES = TRUE)
names(abstracts) <- names(abstracts)  %>% sub(".*abstract_","",.) %>% sub(".csv","",.)
bioprojects[match(names(abstracts),id), abstract := abstracts]

fwrite(bioprojects, "bioprojects.csv")
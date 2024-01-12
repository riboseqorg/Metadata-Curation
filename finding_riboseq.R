library(data.table)
library(ORFik)

studies <- sapply(system("ls SraRunInfo/SraRunInfo*", intern = TRUE), fread, USE.NAMES = TRUE, simplify = FALSE)
studies <- lapply(studies, function(x) as.data.table(apply(x,2, function(y) as.character(y))))
names(studies) <- names(studies) %>%  sub(".*SraRunInfo_","",.) %>% sub(".csv","",.)
dt <- rbindlist(studies, fill = TRUE)
dt <- dt[!duplicated(Run),]

termindexer <-  function(terms, dt) {
  out <- lapply(terms,  function(z) apply(dt,1, function(x) grep(z,x, ignore.case = T)))
  out1 <- lapply(out, function(x) which(0 != (x %>% lapply(length) %>% unlist)))
  out2 <- unlist(out1) %>% unique
  dt[out2]
}

basic_terms <- c("ribo", "footprint", "RPF", "RFP", "80S", "MNAse", "translatome")

first_batch <- termindexer(basic_terms, dt)
second_batch <- termindexer(c("rp", "rf", "fp"),  dt[!(Run %in% first_batch$Run),][,c('Run','LibraryName','sample_title'), ])

# If the samples in the samples whitelist are not in the termindex with the high quality terms then 
# search the second tier list. If found include. If not found then search the main df. 
#If missing from there then there is an issue with the BioProjects whitelist. 
#All samples on the whitelist should have their BioProjects in the in the BioProjects whitelist
fwrite(first_batch, "temp_files/RiboSeq_Metadata_All_Columns.csv")
fwrite(second_batch, "temp_files/RiboSeq_Metadata_second_batch.csv")

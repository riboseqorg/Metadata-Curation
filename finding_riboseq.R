library(data.table)
library(ORFik)

studies <- sapply(system("ls /data/SraRunInfo/SraRunInfo*", intern = TRUE), fread, USE.NAMES = TRUE, simplify = FALSE)
studies <- lapply(studies, function(x) as.data.table(apply(x,2, function(y) as.character(y))))
names(studies) <- names(studies) %>%  sub(".*SraRunInfo_","",.) %>% sub(".csv","",.)

print("Checking studies content:")
print(length(studies))

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
# print(first_batch)
second_batch <- termindexer(c("rp", "rf", "fp"),  dt[!(Run %in% first_batch$Run),][,c('Run','LibraryName','sample_title'), ])

# Path to whitelist CSV - you can change this path as needed
whitelisted_path <- "/usr/local/share/riboseq/resources/whitelisted_samples.csv"

samples_whitelist <- fread(whitelisted_path)

# Identify samples in the whitelist that are not in the first_batch
missing_samples <- setdiff(samples_whitelist$Run, first_batch$Run)
additional_samples <- dt[Run %in% missing_samples]
first_batch <- rbind(first_batch, additional_samples)

fwrite(first_batch, "/data/temp_files/RiboSeq_Metadata_All_Columns.csv")
fwrite(second_batch, "/data/temp_files/RiboSeq_Metadata_second_batch.csv")

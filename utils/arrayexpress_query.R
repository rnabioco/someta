library(tidyverse)
library(ArrayExpress)
library(here)
library(data.table)

# functions
get_samplemeta <- function(id) {
  samplemeta <- paste0("https://www.ebi.ac.uk/arrayexpress/files/",
                       id,
                       "/",
                       id,
                       ".sdrf.txt")
  samplemeta
}

get_proc_files <- function(samplemeta) {
  temp <- fread(samplemeta)
  tech <- temp$`Comment[library construction]`[1]
  files <- as.data.frame(temp)[, which(colnames(temp) == "Derived Array Data File")] %>%
    unlist() %>%
    unique()
  files <- files[!(files == "")] %>% str_c(collapse = ";")
  list(files, tech)
}

get_expmeta <- function(id) {
  expmeta <- paste0("https://www.ebi.ac.uk/arrayexpress/files/",
                    id,
                    "/",
                    id,
                    ".idf.txt")
  expmeta
}

get_add_files <- function(expmeta) {
  temp <- suppressWarnings(readLines(url(expmeta)))
  lines <- temp[str_detect(temp, "\\[AdditionalFile")]
  files <- lines %>% str_remove("..+\\]\t") %>% str_c(collapse = ";")
  files
}

# read query results
message("parse arrayexpress 10x results")

# according to guidelines
sets <- queryAE(keywords = "RNA-seq of coding RNA from single cells")

sets <- sets %>%
  mutate(samplemeta = map_chr(ID, get_samplemeta)) %>%
  mutate(samplemetainfo = map(samplemeta, get_proc_files)) %>%
  unnest_wider(samplemetainfo, names_sep = "_")

sets <- sets %>%
  mutate(expmeta = map_chr(ID, get_expmeta)) %>%
  mutate(additionalfiles = map_chr(expmeta, get_add_files))

write_tsv(sets, here("inst", "extdata", "arrayexpress_sets_proc_add.tsv.gz"))

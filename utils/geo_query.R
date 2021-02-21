library(tidyverse)
library(here)
library(GEOquery)
library(easyPubMed)
library(rcrossref)
date <- format(format = "%m%d%y", Sys.Date())
proj_dir <- here()
fileloc <- file.path(proj_dir, "inst", "extdata", date, paste0("gds_result_", date, ".txt"))

# read query results
message("parse GEO results")
gds <- read_table(fileloc, col_names = FALSE) %>%
  mutate(num = as.numeric(str_remove(str_extract(X1, "^[0-9]+\\."), "\\."))) %>%
  fill(num, .direction = "down") %>%
  filter(str_detect(X1, "Organism:|FTP download:")) %>%
  mutate(type = ifelse(str_detect(X1, "Organism:"), "org", "files")) %>%
  pivot_wider(names_from = type, values_from = X1) %>%
  mutate(id = str_remove(str_remove(files, ".+nnn/"), "/")) %>%
  mutate(org = str_remove(org, "Organism:	")) %>%
  mutate(link = str_c(str_remove(files, ".+\\) "), "suppl/filelist.txt")) %>%
  mutate(link = str_remove(link, "^FTP download: GEO ")) %>%
  mutate(files = str_remove(str_remove(files, "^FTP download: GEO \\("), "\\)..+")) %>%
  mutate(files = ifelse(str_detect(files, "FTP download: GEO"), "none", files)) %>%
  select(id, org, files, link)

# functions
geo_supp <- function(id, string) {
  tryCatch(suppressMessages(getGEOSuppFiles(id,
    makeDirectory = F,
    fetch_files = F
  )),
  error = function(e) {
    "error_get"
  }
  )
}

supp_files <- function(temp) {
  tryCatch(
    temp$fname %>%
      str_to_lower(),
    error = function(e) {
      "none"
    }
  )
}

tar_files <- function(link) {
  tryCatch(
    read_tsv(link, comment = "#") %>% pull(2) %>%
      str_to_lower(),
    error = function(e) {
      "error_parse"
    }
  )
}

check_string <- function(files, string) {
  ifelse(str_detect(files, string) %>% any(),
    "yes",
    "no"
  )
}

get_geo <- function(geo) {
  res <- tryCatch(
    suppressMessages(GEOquery::getGEO(
      GEO = geo,
      filename = NULL,
      GSElimits = NULL, GSEMatrix = FALSE,
      AnnotGPL = FALSE, getGPL = FALSE,
      parseCharacteristics = FALSE
    )),
    error = function(e) {
      "notavail"
    }
  )
  if (class(res) != "GSE") {
    return(NA)
  }
  res@header
}

get_date <- function(res) {
  if (!is.na(res)) {
    res$submission_date
  } else {
    NA
  }
}

get_pubmed <- function(res) {
  # res <- res[[1]]
  if (class(res) != "list") {
    return(NA)
  }

  if (length(res$pubmed_id) > 1) {
    id2 <- tryCatch(res$pubmed_id[1],
      error = function(e) {
        "notfound"
      }
    )
  } else {
    id2 <- tryCatch(res$pubmed_id,
      error = function(e) {
        "notfound"
      }
    )
  }

  if (is.null(id2)) {
    return(NA)
  }
  if (id2 == "notfound") {
    return(NA)
  }
  return(id2)
  # tryCatch(
  #   get_pubmed_ids(id2) %>%
  #     fetch_pubmed_data() %>%
  #     article_to_df(),
  #   error = function(e) {"error"}
  # )
}

from_list <- function(pubmedid) {
  tryCatch(pubmeds[[pubmedid]],
    error = function(e) {
      "error"
    }
  )
}

get_journal <- function(pubmed) {
  tryCatch(
    pubmed %>%
      dplyr::slice(1) %>%
      pull(jabbrv),
    error = function(e) "error"
  )
}

# step-by-step processing
message("parse GEO supplementary files")
gds <- gds %>% mutate(supp = pbmcapply::pbmcmapply(
  geo_supp,
  id
))
gds2 <- gds %>% mutate(suppfiles = pbmcapply::pbmcmapply(
  supp_files,
  supp
))
gds3 <- gds2 %>% mutate(tarfiles = pbmcapply::pbmcmapply(
  tar_files,
  link
))
gds4 <- gds3 %>% mutate(has_meta1 = pbmcapply::pbmcmapply(
  check_string,
  suppfiles, "meta|annot|type|clustering|coldata"
))
gds4 <- gds4 %>% mutate(has_meta2 = pbmcapply::pbmcmapply(
  check_string,
  tarfiles, "meta|annot|type|clustering|coldata"
))
gds4 <- gds4 %>% mutate(has_object1 = ifelse(str_detect(suppfiles, "rds|rda|rdata|loom|h5ad"),
  "yes", "no"
))
gds4 <- gds4 %>% mutate(has_object2 = ifelse(str_detect(tarfiles, "rds|rda|rdata|loom|h5ad"),
  "yes", "no"
))
gds4 <- gds4 %>% mutate(usable = ifelse(has_meta1 == "yes" | has_meta2 == "yes" | has_object1 == "yes" | has_object2 == "yes",
  "yes", "no"
))

# giant superseries (also incorrect) soft.gz  files cause issues
blacklist <- c("GSE63058", "GSE47917")
gds4 <- gds4 %>% filter(!(id %in% blacklist))

# geo and pubmed parsing, slow
message("get GEO data")
gds5 <- gds4 %>%
  split(gds4$id) %>%
  map_dfr(. %>% mutate(geo = map(id, get_geo)))
saveRDS(gds5, here("inst", "extdata", date, paste0("gds5_temp_", date, ".rds")))

message("get PubMed data")
gds5 <- readRDS(here("inst", "extdata", date, paste0("gds5_temp_", date, ".rds")))
gds5 <- gds5 %>%
  mutate(date = pbmcapply::pbmcmapply(get_date, geo)) %>%
  separate(date, into = c("month", "day", "year"))
gds5 <- gds5 %>% mutate(pubmed = pbmcapply::pbmcmapply(get_pubmed, geo))

pubmeds <- list()
splits <- split(gds5$pubmed %>% na.omit(), ceiling(seq_along(gds5$pubmed %>% na.omit()) / 20))
for (s in splits) {
  query_string <- str_c(s, collapse = ",")
  pubmed_file <- get_pubmed_ids(query_string,
    api_key = "9cc836a53f5c509a6783283e694da4fb0d08"
  )
  pubmed <- pubmed_file %>%
    fetch_pubmed_data(retmax = 9999) %>%
    articles_to_list()
  pubmeds1 <- sapply(pubmed, article_to_df, simplify = FALSE)
  names(pubmeds1) <- pubmed_file$IdList
  pubmeds <- append(pubmeds, pubmeds1)
}

gds6 <- gds5 %>% mutate(pubmed = map(pubmed, from_list))
gds6 <- gds6 %>% mutate(journal = pbmcapply::pbmcmapply(
  get_journal,
  pubmed
))

get_cites <- function(pubmed) {
  tryCatch(
    cr_citation_count(pubmed$doi[1])$count,
    error = function(e) "NA"
  )
}
message("get citation data")
gds7 <- gds6 %>%
  mutate(cite = map(pubmed, get_cites)) %>%
  mutate(year = as.numeric(year))
gds7 %>% saveRDS(here("inst", "extdata", date, paste0("geo_", date, ".rds")))
gds7 %>% saveRDS(here("inst", "extdata", "current_geo.rds"))

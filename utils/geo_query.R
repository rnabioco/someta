library(tidyverse)
library(here)
library(GEOquery)
proj_dir <- here()
# this file of exported geo search results
# "expression profiling by high throughput sequencing"[DataSet Type] AND "single cell"[All Fields] 
fileloc <- file.path(proj_dir, "inst", "extdata", paste0("gds_result_", date, ".txt"))

gds <- read_table(fileloc, col_names = FALSE) %>% 
  mutate(num = as.numeric(str_remove(str_extract(X1, "^[0-9]+\\."), "\\."))) %>% 
  fill(num, .direction = "down") %>% 
  filter(str_detect(X1, "Organism:|FTP download:")) %>% 
  mutate(type = ifelse(str_detect(X1, "Organism:"), "org", "files")) %>% 
  pivot_wider(names_from = type, values_from = X1) %>% 
  mutate(id = str_remove(str_remove(files, ".+nnn/"), "/")) %>% 
  mutate(org = str_remove(org, "Organism:	")) %>% 
  mutate(files = str_remove(str_remove(files, "^FTP download: GEO \\("), "\\)..+")) %>% 
  mutate(files = ifelse(str_detect(files, "FTP download: GEO"), "none", files)) %>% 
  select(id, org, files)

geo_string <- function(id, string) {
  temp <- tryCatch(suppressMessages(getGEOSuppFiles(id, 
                                                    makeDirectory = F, 
                                                    fetch_files = F)),
                   error = function(e) {"error_get"})
  if (is.null(temp)) {
    return("no")
  } else {
    tryCatch(
      ifelse(
        temp$fname %>% 
          str_to_lower() %>%
          str_detect(string) %>% 
          any(),
        "yes",
        "no"),
      error = function(e) {"error_parse"})
  }
}

t0 <- Sys.time()
gds <- gds %>%
  mutate(has_meta = pbmcapply::pbmcmapply(geo_string, id, "meta")) %>%
  filter(!str_detect(has_meta, "error")) %>%
  mutate(has_r = ifelse(str_detect(files, "RDA|RDATA|RDS"), "yes", "no")) %>% 
  mutate(usable = ifelse(has_meta == "yes" | has_r == "yes", "yes", "no"))
message("GEOquery data check step took ", format(Sys.time() - t0))

gds_h <- gds %>% filter(str_detect(org, "Homo sapiens"))
gds_m <- gds %>% filter(str_detect(org, "Mus musculus"))

write_tsv(gds_h, here("inst", "extdata", paste0("geo_hs_", date, ".tsv")))
write_tsv(gds_m, here("inst", "extdata", paste0("geo_mm_", date, ".tsv")))
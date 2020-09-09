library(tidyverse)
library(here)
proj_dir <- here()
# this file of exported geo search results
# "expression profiling by high throughput sequencing"[DataSet Type] AND "single cell"[All Fields] 
fileloc <- file.path(proj_dir, "data", "geo_query", "gds_result_allsc.txt")

gds <- read_table(fileloc, col_names = FALSE) %>% 
  filter(str_detect(X1, "Organism:|FTP download:")) %>% 
  mutate(entry = floor((row_number() - 1) / 2), part = (row_number() - 1) %% 2) %>% 
  pivot_wider(names_from = part, values_from = X1) %>% 
  mutate(id = str_remove(str_remove(`1`, ".+nnn/"), "/")) %>% 
  mutate(org = str_remove(`0`, "Organism:	")) %>% 
  mutate(files = str_remove(str_remove(`1`, "^FTP download: GEO \\("), "\\)..+")) %>% 
  filter(!str_detect(files, "FTP download: GEO")) %>% 
  select(id, org, files)

geo_string <- function(id, string) {
  temp <- tryCatch(suppressMessages(getGEOSuppFiles(id, 
                                                    makeDirectory = F, 
                                                    fetch_files = F)),
                   error = function(e) {"error1"})
  if (is.null(temp)) {
    return(FALSE)
  } else {
    tryCatch(
      temp$fname %>% 
        str_to_lower() %>%
        str_detect(string) %>% 
        any(),
      error = function(e) {"error2"})
  }
}

gds_h <- gds %>% filter(str_detect(org, "Homo sapiens")) %>%
  mutate(has_meta = pbmcapply::pbmcmapply(geo_string, id, "meta"))
gds_h <- gds_h %>% mutate(has_r = str_detect(files, "RDA|RDATA|RDS")) %>% 
  mutate(usable = ifelse(has_meta | has_r, TRUE, FALSE))

gds_m <- gds %>% filter(str_detect(org, "Mus musculus")) %>%
  mutate(has_meta = pbmcapply::pbmcmapply(geo_string, id, "meta"))
gds_m <- gds_m %>% filter(!str_detect(has_meta, "error")) %>%
  mutate(has_meta = as.logical(has_meta)) %>% 
  mutate(has_r = str_detect(files, "RDA|RDATA|RDS")) %>% 
  mutate(usable = ifelse(has_meta | has_r, TRUE, FALSE))

write_tsv(gds_h, "geo_hs.tsv")
write_tsv(gds_m, "geo_mm.tsv")
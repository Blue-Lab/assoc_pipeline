#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Add FID to sample ID") %>%
  add_argument("plink_prefix", help  = "Prefix for PLINK bed/bim/fam files") %>%
  add_argument("sample_id", help = "Vector of sample IDs (.rds)") %>%
  add_argument("--out_file", default = "sample_ids.txt", help = "output filename")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(dplyr)

sample.id <- readRDS(argv$sample_id)
fam_header <- c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHENOTYPE")
fam <- paste0(argv$plink_prefix, ".fam") %>% read.table(col.names = fam_header)

fam_keep <- filter(fam, IID %in% sample.id) %>% select(FID, IID)
write.table(fam_keep, argv$out_file, na = "", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

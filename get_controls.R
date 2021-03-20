#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Get IDs of controls") %>%
  add_argument("pheno_file", help="Phenotype file in annotated dataframe format (.rds)") %>%
  add_argument("status", help = "Case/control status variable in pheno_file (binary)") %>%
  add_argument("plink_prefix", help  = "Prefix for PLINK bed/bim/fam files") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs to keep (.rds)")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(Biobase)
library(dplyr)

fam_header <- c("FID", "IID", "FATHER", "MOTHER", "SEX", "PHENOTYPE")
fam <- paste0(argv$plink_prefix, ".fam") %>% read.table(col.names = fam_header)

pheno <- readRDS(argv$pheno_file) %>% pData()
controls <- filter(pheno, .data[[argv$status]] == 0)$sample.id

if (!is.na(argv$sample_id)) {
  sample_id <- readRDS(argv$sample_id)
  keep <- controls[controls %in% sample_id]
} else {
  keep <- controls
}

fam_keep <- filter(fam, IID %in% keep) %>% select(FID, IID)
paste0(argv$out_prefix, "controls.txt") %>%
  write.table(fam_keep, ., na = "", row.names = FALSE, col.names = FALSE,
              quote = FALSE)

paste0(argv$out_prefix, "controls.rds") %>% saveRDS(keep, .)

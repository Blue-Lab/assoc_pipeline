#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Format PC-AiR PCs for RUTH") %>%
  add_argument("pc_file", help = "File with pcair PCs object (.rds)") %>%
  add_argument("fam_file", help = "PLINK .fam file. Can be partial version from get_controls.R") %>%
  add_argument("--out_file", help = "Name for output file",
               default = "ruth_pcs.txt")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(dplyr)

fam <- read.table(argv$fam_file)[,1:2] %>%
  set_colnames(c("FID", "SAMPLE_ID"))

pcair <- readRDS(argv$pc_file)
pcs <- readRDS(argv$pc_file) %>%
  as.data.frame() %>%
  # variable names must be SAMPLE_ID, PC1, PC2, ....
  rename_with(~ gsub("V", "PC", .x)) %>%
  mutate(SAMPLE_ID = rownames(.)) %>%
  left_join(fam, by = "SAMPLE_ID") %>%
  # RUTH requires the ID field to be a concatenation of FID and IID
  mutate(SAMPLE_ID = paste(FID, SAMPLE_ID, sep = "_")) %>%
  select(SAMPLE_ID, matches("PC"))
write.table(pcs, argv$out_file, sep = "\t", na = "", row.names = FALSE,
            quote = FALSE)

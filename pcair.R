#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Run PC-AiR") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("kin_file", help = "Kinship matrix (KING or PC-Relate)") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--out_file", help = "output file name",
               default = "pcair.rds") %>%
  add_argument("--kin_thresh", help = "Kinship threshold for pcair",
               default = 0)

argv <- parse_args(argp)

library(SeqArray)
library(GENESIS)

sessionInfo()
print(argv)

gds <- seqOpen(argv$gds_file)

if (!is.na(argv$variant_id)) {
  variant_id <- readRDS(argv$variant_id)
} else {
  variant_id <- NULL
}
if (!is.na(argv$variant_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}

kin <- readRDS(argv$kin_file)

if (class(kin) == "snpgdsIBDClass") {
  mat <- kin$kinship
} 

mypcair <- pcair(gds, kinobj = mat, kin.thresh = argv$kin_thresh,
                 divobj = mat, snp.include = variant_id,
                 sample.include = sample_id)

saveRDS(mypcair, argv$out_file)

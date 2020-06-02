#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Generate KING GRM") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("--out_file", help = "output file name",
               default = "king_grm.rds") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(SNPRelate)
library(SeqArray)

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

gds <- seqOpen(argv$gds_file)

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id)

rownames(king$kinship) <- king$sample.id
colnames(king$kinship) <- king$sample.id

saveRDS(king, argv$out_file)

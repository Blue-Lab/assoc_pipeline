#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Run KING-robust") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--num_core", help = "num.thread argument for snpgdsIBDKING (if NA, detect the number of cores automatically)")
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

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id,
                      type = "KING-robust",
                      num.thread = as.numeric(argv$num_core))

rownames(king$kinship) <- king$sample.id
colnames(king$kinship) <- king$sample.id

saveRDS(king, paste0(argv$out_prefix, "king_out.rds"))
saveRDS(king$kinship, paste0(argv$out_prefix, "king_grm.rds"))

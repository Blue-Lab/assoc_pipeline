#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Run KING-robust") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--family_id", help = "File with vector of family IDs") %>%
  add_argument("--num_core", help = "num.thread argument for snpgdsIBDKING (if NA, detect the number of cores automatically)")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(SNPRelate)
library(SeqArray)
library(dplyr)

if (!is.na(argv$variant_id)) {
  variant_id <- readRDS(argv$variant_id)
} else {
  variant_id <- NULL
}

if (!is.na(argv$sample_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}

gds <- seqOpen(argv$gds_file)
gds.sample.id <- seqGetData(gds, "sample.id")

if (!is.na(argv$family_id)) {
  family_id <- readRDS(argv$family_id)
  if (is.null(sample_id)) {
    stop("No --sample_id supplied. Supply --sample_id when using --family_id")
  } else {
    tmp_fid <- data.frame(FID = family_id, sample.id = sample_id) %>%
      filter(sample.id %in% gds.sample.id)
    family_id <- tmp_fid$FID
    sample_id <- tmp_fid$sample.id
  }
} else {
  family_id <- NULL
}

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id,
                      family.id = family_id, type = "KING-robust",
                      num.thread = as.numeric(argv$num_core))

rownames(king$kinship) <- king$sample.id
colnames(king$kinship) <- king$sample.id

saveRDS(king, paste0(argv$out_prefix, "king_out.rds"))
saveRDS(king$kinship, paste0(argv$out_prefix, "king_grm.rds"))

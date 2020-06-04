#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Run PC-Relate") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pcair_file", help = "File with pcair object") %>%
  add_argument("--out_file", help = "output file name",
               default = "pcrelate.rds") %>%
  add_argument("--n_pcs", default = 10,
               "Number of PCs to pass to pcrelate") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(SeqArray)
library(SeqVarTools)
library(GENESIS)

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
mypcair <- readRDS(argv$pcair_file)

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
seqData <- SeqVarData(gds)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(argv$n_pcs)],
                    training.set = mypcair$unrels, sample.include = sample_id)

saveRDS(mypcrel, argv$out_file)

#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

st_help <- "Threshold for sparsifying GRM (will be multiplied by scale_kin)"
argp <- arg_parser("Run PC-Relate") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pcair_file", help = "File with pcair object") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--n_pcs", default = 10,
               "Number of PCs to pass to pcrelate") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--scale_kin", help = "Scaling factor for GRM output",
               default = 1) %>%
  add_argument("--sparse_thresh", help = "st_help", default = 0) %>%
  add_argument("--scan_block_size", help = "pcrelate scan block size",
               default = 5000) %>%
  add_argument("--correct", "Flag to implement small-sample correction",
               flag = TRUE) %>%
  add_argument("--variant_block", help = "SeqVarBlaockIterator block size",
               default = 1024)
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
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE,
                                variantBlock = argv$variant_block)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(argv$n_pcs)],
                    training.set = mypcair$unrels, sample.include = sample_id,
                    scan.block.size = argv$scan_block_size)

saveRDS(mypcrel, paste0(argv$out_prefix, "pcrelate.rds"))
pcr_mat <- pcrelateToMatrix(mypcrel, thresh = argv$sparse_thresh,
                            scaleKin = argv$scale_kin, correct = argv$correct)
saveRDS(pcr_mat, paste0(argv$out_prefix, "pcr_mat.rds"))

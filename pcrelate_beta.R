#! /usr/bin/env Rscript

library(argparser)
library(magrittr)

argp <- arg_parser("PC-Relate") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pca_file", help = "PC-AiR output file (.rds)") %>%
  add_argument("n_pcs", help = "Number of PCs to include") %>%
  add_argument("--out_file", help = "Name for output file", default = "pcrelate_beta.rds") %>%
  add_argument("--sample_id", help = "Vector of sample IDs to include (.rds)") %>%
  add_argument("--variant_id", help = "Vector of variant IDs to include (.rds)") %>%
  add_argument("--variant_block_size", default = 1024)
argv <- parse_args(argp)

library(SeqVarTools)
library(GENESIS)

sessionInfo()
print(argv)

gds <- seqOpen(argv$gds_file)
seqData <- SeqVarData(gds)

if (!is.na(argv$variant_include_file)) {
  variant_id <- readRDS(argv$variant_id)
  seqSetFilter(gds, variant.id = variant_id, action = "intersect")
}

pca <- readRDS(argv$pca_file)
n_pcs <- min(as.integer(argv$n_pcs), length(pca$unrels))
pcs <- as.matrix(pca$vectors[,1:n_pcs])
sample.include <- samplesGdsOrder(seqData, pca$unrels)

if (!is.na(argv$sample_id)) {
  sample.id <- readRDS(argv$sample_id)
  sample.include %<>% intersect(sample.id)
}


# create iterator
block.size <- as.integer(argv$variant_block_size)
iterator <- SeqVarBlockIterator(seqData, variantBlock = block.size)

beta <- calcISAFBeta(iterator, pcs = pcs, sample.include = sample.include)

saveRDS(beta, argv$out_file)

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")

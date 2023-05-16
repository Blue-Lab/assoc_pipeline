#! /usr/bin/env Rscript

library(argparser)
library(magrittr)

argp <- arg_parser("PC-Relate by bloc") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pca_file", help = "PC-AiR output (.rds)") %>%
  add_argument("beta_file", help = "ISAF beta file (.rds)") %>%
  add_argument("block_1", help = "Sample-include block 1", type = "integer") %>%
  add_argument("block_2", help = "Sample-include block 2", type = "integer") %>%
  add_argument("n_sample_blocks", help = "Number of sample blocks", type = "integer") %>%
  add_argument("n_pcs", help = "Number of PCs to include", type = "integer") %>%
  add_argument("--out_prefix", help = "Prefix for output filename", default = "pcrelate") %>%
  add_argument("--sample_id", help = "Vector of sample IDs (.rds)") %>%
  add_argument("--variant_block_size", help = "Block size for SeqVarBlockIterator", default = 1024L) %>%
  add_argument("--variant_id", help = "Vector of variant IDs (.rds)") %>%
  add_argument("--ibd_probs", help = "Estimate pairwise IBD sharing probabilities?", default = TRUE)
argv <- parse_args(argp)

library(SeqVarTools)
library(GENESIS)

sessionInfo()
print(argv)
print(dput(argv))
n_sample_blocks <- as.integer(argv$n_sample_blocks)
i <- as.integer(argv$block_1)
j <- as.integer(argv$block_2)

maxblock <- max(i, j)
if (maxblock > n_sample_blocks) {
  "You have specified there are only %s total blocks, but passed block #%s" %>%
    sprintf(n_sample_blocks, maxblock) %>%
    stop()
}

gds <- seqOpen(argv$gds_file)
seqData <- SeqVarData(gds)
if (!is.na(argv$variant_id)) seqSetFilter(seqData, readRDS(argv$variant_id))

pca <- readRDS(argv$pca_file)
n_pcs <- min(as.integer(argv$n_pcs), length(pca$unrels))
pcs <- as.matrix(pca$vectors[,1:n_pcs])
sample.include <- samplesGdsOrder(seqData, rownames(pcs))

if (!is.na(argv$sample_id)) {
    sample.id <- readRDS(argv$sample_id)
    sample.include <- intersect(sample.include, sample.id)
}

# load betas
betaobj <- readRDS(argv$beta_file)

# create iterator
iterator <- SeqVarBlockIterator(seqData, variantBlock = argv$variant_block_size)

# create sample blocks
if (n_sample_blocks > 1) {
    samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), n_sample_blocks)))
} else {
    samp.blocks <- list(sample.include)
}

out <- pcrelateSampBlock(iterator, betaobj = betaobj, pcs = pcs,
                         sample.include.block1 = samp.blocks[[i]],
                         sample.include.block2 = samp.blocks[[j]],
                         ibd.probs = argv$ibd_probs)

saveRDS(out, file = paste0(argv$out_prefix, "_block_", i, "_", j, ".rds"))

seqClose(seqData)

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")

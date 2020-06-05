#! /usr/bin/env Rscript
library(argparser)
msg <- "This script is defunct. Consider using the scripts king_grm.R," %>%
  paste("pcair.R and pcrelate.R")
warning(msg)
library(magrittr)

# read arguments
argp <- arg_parser("Generate PCs and GRM") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--kin_thresh", default = 5.5,
               help = "Kinship threshold for pcair (2 ^ -kin_thresh)") %>%
  add_argument("--n_pcs", default = 3,
               "Number of PCs to pass to pcrelate") %>%
  add_argument("--keep_king", flag = TRUE, help = "Save KING-robust GRM")
argv <- parse_args(argp)

library(SeqArray)
library(GENESIS)
library(SeqVarTools)
library(SNPRelate)

sessionInfo()
print(argv)

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
kin_thresh <- 2 ^ (-argv$kin_thresh)
out_prefix <- argv$out_prefix
gds <- seqOpen(argv$gds_file)

king <- snpgdsIBDKING(gds, snp.id = variant_id, sample.id = sample_id)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

if (argv$keep_king) {
  kingMat_temp <- kingMat * 2 # Scaled to match pc-relate GRM
  # coerces low values in matrix to 0
  kingMat_temp[kingMat_temp <= kin_thresh] <- 0
  saveRDS(kingMat_temp, paste0(out_prefix, "king_robust_grm.rds"))
  rm(kingMat_temp)
}

mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = kin_thresh,
                 divobj = kingMat, snp.include = variant_id,
                 sample.include = sample_id)

seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
seqData <- SeqVarData(gds)
print("1st iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, seq(argv$n_pcs)],
                    training.set = mypcair$unrels)
pcrelate_matrix <- pcrelateToMatrix(mypcrel, scaleKin=2, thresh = kin_thresh)

pca <- pcair(seqData, kinobj = pcrelate_matrix, kin.thresh = kin_thresh,
             divobj = kingMat, snp.include = variant_id,
             sample.include = sample_id)

resetIterator(iterator, verbose = TRUE)

print("2nd iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)
pcrel2 <- pcrelate(iterator, pcs = pca$vectors[, seq(argv$n_pcs)],
                   training.set = pca$unrels)

pcrelate_matrix2 <- pcrelateToMatrix(pcrel2, scaleKin = 2, thresh = kin_thresh)
saveRDS(pca, paste0(out_prefix, "pcair.rds"))
saveRDS(pcrelate_matrix2, paste0(out_prefix, "pcr_grm.rds"))

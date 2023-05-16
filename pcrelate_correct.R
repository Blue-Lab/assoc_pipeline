#! /usr/bin/env Rscript

library(argparser)
library(magrittr)

argp <- arg_parser("PC-Relate: combine and correct blocks") %>%
  add_argument("pcrelate_prefix", help = "prefix for pcrelate_block.R output files") %>%
  add_argument("n_sample_blocks",
               help = "Number of blocks to divide samples into") %>%
  add_argument("--sparse_thresh", default = 0,
               help = "Threshold for sparsifying GRM (will be multiplied by scale_kin)") %>%
  add_argument("--scale_kin", help = "Scaling factor for GRM output",
               default = 1)

argv <- parse_args(argp)
library(GENESIS)
library(data.table)

sessionInfo()
print(argv)

nsampblock <- as.integer(argv$n_sample_blocks)

kinSelf <- NULL
kinBtwn <- NULL
kin.thresh <- as.numeric(argv$sparse_thresh)

# correct IBD results and combine
# Note the index in the outer loop is reversed from the usual
for (i in nsampblock:1){
    for (j in i:nsampblock){
        message('Sample Blocks ', i, ' and ', j)
        
        ## load the data
        res <- readRDS(paste0(argv$pcrelate_prefix, "_block_", i, "_", j, ".rds")) 
        
        if (i == j) kinSelf <- rbind(kinSelf, res$kinSelf)
        
        # correct the IBD estimates
        res$kinBtwn <- correctK2(kinBtwn = res$kinBtwn, 
                                 kinSelf = kinSelf, 
                                 small.samp.correct = FALSE, 
                                 pcs = NULL, 
                                 sample.include = NULL)
        
        res$kinBtwn <- correctK0(kinBtwn = res$kinBtwn)
        saveRDS(res, file = paste0(argv$pcrelate_prefix, "_block_", i, "_", j, "_corrected.rds"))

        # save results above threshold in combined file
        kinBtwn <- rbind(kinBtwn, res$kinBtwn[kin > kin.thresh])
        
        rm(res); gc()
    }
}

# save pcrelate object
pcrelobj <- list(kinSelf = kinSelf, kinBtwn = kinBtwn)
class(pcrelobj) <- "pcrelate"
saveRDS(pcrelobj, file = paste0(argv$pcrelate_prefix, "_pcrelate.rds"))

rm(kinBtwn, kinSelf); gc()
   
# save sparse kinship matrix
km <- pcrelateToMatrix(pcrelobj, thresh = argv$scale_kin * kin.thresh, scaleKin = argv$scale_kin)
saveRDS(km, file = paste0(argv$pcrelate_prefix, "_pcrelate_Matrix.rds"))

# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")

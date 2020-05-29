#! /usr/bin/env Rscript
#' Calculate lambda, the genomic inflation factor
#'
#' @param stat Vector of test statistics
#' @param df Degrees of freedom
#' @return lambda
#'
#' @importFrom stats median qchisq
#' @export
calculateLambda <- function(stat, df) {
    if (any(sum(stat < 0, na.rm=TRUE)))
        stop("no negative values allowed in stat (does beta/se need to be squared?)")
    median(stat, na.rm=TRUE) / qchisq(0.5, df=df)
}

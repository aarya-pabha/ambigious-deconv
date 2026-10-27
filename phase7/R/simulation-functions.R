#' Makes a binary endmember matrix
#' @param m The number of snps per sequence
#' @param r The number of sequences
#' @param full_rank If TRUE, will create an endmember matrix of full
#'     rank. If FALSE, will create an endmember matrix of rank
#'     r-1. full_rank = FALSE only allowed if r >= 3.
#' @return An m x r binary matrix.
make_endmembers <- function(m, r, full_rank = TRUE) {
    W_true = matrix(sample(c(0,1), size = m * r, replace = TRUE), nrow = m, ncol = r)
    all_zero = which(rowSums(W_true) == 0)
    all_one = which(rowSums(W_true) == r)
    for(idx in all_zero) {
        W_true[idx, sample(1:r, size = 1)] = 1
    }
    for(idx in all_one) {
        W_true[idx, sample(1:r, size = 1)] = 0
    }
    W_true
}

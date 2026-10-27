## source("nmf-functions.R")
## ## the following doesn't get the right sequences, but larger values of N do.
## set.seed(22)
## m <- 5; n <- 5; r <- 3; N <- 4


## ## True W (binary endmembers)
## (W_true <- matrix(sample(c(0,1), size = m * r, replace = TRUE), nrow = m, ncol = r))
## all_zero = which(rowSums(W_true) == 0)
## all_one = which(rowSums(W_true) == r)
## for(idx in all_zero) {
##     W_true[idx, sample(1:r, size = 1)] = 1
## }
## for(idx in all_one) {
##     W_true[idx, sample(1:r, size = 1)] = 0
## }
## W_true

## # Dirichlet H and proportions S
## H_true <- matrix(0, r, n)
## for (j in 1:n) H_true[, j] <- as.numeric(gtools::rdirichlet(1, rep(1, r)))
## S <- W_true %*% H_true
## # Observations: binomial counts with per-column trials N_j
## X <- matrix(rbinom(m*n, size = N, prob = as.vector(S)), m, n)

## rays_out <- build_candidates_rays(X, N, r + 1,
##                                   transform = "arcsin",
##                                   n_rays = 1000)
## rays_out$candidates

## W_true_str = apply(W_true, 2, paste0, collapse = "")
## candidates_str = apply(rays_out$candidates, 2, paste0, collapse = "")
## sapply(W_true_str, function(x) x %in% candidates_str)


## out = truncated_posterior_W(X, N = rep(N, ncol(X)), r = r , Wcat = rays_out$candidates)

## out$posterior |> head()
## canon_key(W_true) %in% out$posterior$key
## out$posterior[out$posterior$key == canon_key(W_true),]

## Wmap = out$posterior$key[1] |> canonical_key_to_W()
## Phat = diag_fit_WH(X, rep(N, ncol(X)), Wmap, eps = 1e-12)$Phat

## Phat2 = out$posterior$key[2] |> canonical_key_to_W() |> diag_fit_WH(X, rep(N, ncol(X)),W =  _) |> extract2("Phat")
## Phat1 = out$posterior$key[1] |> canonical_key_to_W() |> diag_fit_WH(X, rep(N, ncol(X)),W =  _) |> extract2("Phat")
## Phat3 = out$posterior$key[3] |> canonical_key_to_W() |> diag_fit_WH(X, rep(N, ncol(X)),W =  _) |> extract2("Phat")

## boot_out = bootstrap_no_refit_loglik(X, N, Phat, B = 1000, eps = 1e-12)
## boot_out$p_value


## bootstrap_no_refit_loglik(X, N, Phat1, B = 1000, eps = 1e-12)$p_value
## bootstrap_no_refit_loglik(X, N, Phat2, B = 1000, eps = 1e-12)$p_value
## bootstrap_no_refit_loglik(X, N, Phat3, B = 1000, eps = 1e-12)$p_value

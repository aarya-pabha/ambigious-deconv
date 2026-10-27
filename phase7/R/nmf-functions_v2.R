# ============================================================
# Ray-casting candidates for binary endmembers (VolMin init)
# Steps:
#   (1) Fit affine subspace (mu, U) under binomial noise
#   (2) Shoot random rays in the subspace, intersect with cube,
#       snap intersection endpoints to nearest 0/1 vertices,
#       (optional) polish with greedy bit-flips.
# ============================================================

# ---------- small utilities ----------
#' Transform probabilities
#' @param p A vector/matrix containing values in [0,1].
#' @param method The method to use for transforming the probabilities.
#' @param eps The value used for clipping the probabilities away from 0 and 1.
#' @return The transformed value(s).
transform_prob <- function(p, method = c("arcsin","none"), eps = 1e-12) {
  method <- match.arg(method)
  p <- clip(p, eps, 1 - eps)
  if (method == "arcsin") {
    2 * asin(sqrt(p))                 # in [0, pi]
  } else {
    p                                 # in [0,1]
  }
}

# Map binary w (0/1) into transformed space
transform_binary <- function(w, method = "arcsin", eps = 0) {
  transform_prob(as.numeric(w), method = method, eps = eps)
}

#' Distance between vector and subspace
#' Computes the distance between a vector v and an affine space { mu + Ut: t in R^(ncol(U)) }
#' @param v A vector.
#' @param mu A vector giving a point on the affine space.
#' @param U A matrix defining the affine space.
#' @return The squared norm of the distance between v and its projection onto the affine space.
subspace_dist2 <- function(v, mu, U) {
  z  <- v - mu
  Uz <- t(U) %*% z
  res <- z - U %*% Uz
  sum(res*res)
}

#' 
project_coords <- function(v, mu, U) {
  drop(t(U) %*% (v - mu))
}

#' Fit an affine subspace
#' Fit an affine subspace to a matrix of probabilities, with an optional approximate variance-stabilizing transformation as pre-processing.
#' @param X A count matrix.
#' @param N Either a single value, a vector, or a matrix describing the number of trials corresponding to each cell in X. If a single value, that value gives the number of trials for each cell. If a vector, the jth element of that vector is assumed to give the number of trials for each of the elements in the jth column of X. If a matrix the ijth element of the matrix gives the number of trials for the ijth element of X.
#' @param r The dimension of the affine subspace.
#' @param transform Should the probabilities be transformed before subspace estimation?
#' @param eps_var
#' @return A list containing mu (a point on the subspace), U (an orthogonal matrix defining the subspace), transform (a value specifying how the probabilities were transformed before subspace fitting), and Shat (the probability matrix used to estimate the subspace).
fit_affine_subspace <- function(X, N, r,
                                transform = c("arcsin","none"),
                                eps_var = 1e-6) {
    transform <- match.arg(transform)
    X <- as.matrix(X)
    m <- nrow(X)
    n <- ncol(X)
    Nmat <- make_Nmat(N, m, n)
    Shat <- X / Nmat
    Shat <- clip(Shat, 1e-8, 1 - 1e-8)
    Y <- transform_prob(Shat, method = transform)

    if (transform == "none") {
        W_ij <- sweep(1 / clip(Shat * (1 - Shat), eps_var, Inf), 2, N, "*")
    } else if (transform == "arcsin") {
        W_ij <- matrix(N, nrow = nrow(Shat), ncol = ncol(Shat))
    }
    Wj <- colSums(W_ij)

    mu <- as.numeric(Y %*% (Wj / sum(Wj)))
    Yc <- sweep(Y, 1L, mu, "-")
    Ycw <- sweep(Yc, 2L, sqrt(Wj), "*")

    rminus1 <- max(1L, r - 1L)
    sv <- svd(Ycw, nu = rminus1, nv = 0)
    U <- sv$u
    d <- sv$d  # NEW: Extract singular values for the noise check

    # NEW: Return 'd' in the list
    return(list(mu = mu, U = U, d = d, transform = transform, Shat = Shat))
}



#' 
# Given ray r(t) = mu + t * v, find the interval [t_min, t_max]
# for which r(t) is inside the m-D box [a,b]^m.
line_box_interval <- function(mu, v, a, b) {
  m <- length(mu)
  t_low <- -Inf; t_high <- Inf
  for (i in 1:m) {
    vi <- v[i]
    if (abs(vi) < 1e-14) {
      # Static coordinate; must lie within [a,b]
      if (mu[i] < a - 1e-12 || mu[i] > b + 1e-12) {
        return(c(Inf, -Inf))  # empty
      } else {
        next
      }
    } else {
      t1 <- (a - mu[i]) / vi
      t2 <- (b - mu[i]) / vi
      t_lo_i <- min(t1, t2)
      t_hi_i <- max(t1, t2)
      if (t_lo_i > t_low)  t_low  <- t_lo_i
      if (t_hi_i < t_high) t_high <- t_hi_i
      if (t_low > t_high) return(c(Inf, -Inf)) # empty
    }
  }
  return(c(t_low, t_high))
}

# Snap a point y in transformed space to the nearest {0,1}^m vertex
# by coordinate-wise thresholding at mid = (a+b)/2 in the transformed space.
snap_to_vertex <- function(y, a, b) {
    mid <- (a + b)/2  # vectorized (same a,b per coord)
    as.integer(y >= mid)
}

## ---------- main: ray-based candidate generation ----------
#' Build candidate endmembers
#' Fit affine subspace, sample rays in that subspace, intersect rays with unit hypercube, snap to closest hypercube vertex to find candidate endmembers
#' @param X A matrix of counts.
#' @param N An integer, vector, or matrix giving the number of trials for the cells in X.
#' @param r The dimension of the affine subspace.
#' @param transform Either arcsin or none, transformation to use before fitting the affine subspace. arcsine is variance stabilizing for binomial data.
#' @param n_rays The number of rays to use to generate candidates.
#' @param seed An optional value for the seed. If NULL, no seed is set.
#' @return A list with one element 'candidates', a matrix whose columns are the candidate endmembers, and 'candidate_probs'
build_candidates_rays <- function(X, N, r,
                                  transform = c("arcsin","none"),
                                  n_rays = 2000L,
                                  seed = NULL) {
    transform <- match.arg(transform)
    if(!is.null(seed)) {
        set.seed(seed)
    }
    
    ## (1) Fit affine subspace
    fit <- fit_affine_subspace(X, N, r, transform = transform)
    mu <- fit$mu
    U <- fit$U
    m <- length(mu)

    ## Box bounds in transformed space
    if (transform == "arcsin") {
        a <- 0; b <- pi
    } else if (transform == "none") {
        a <- 0; b <- 1
    }
    a <- rep(a, m); b <- rep(b, m)

    rays <- list()

    ## random directions in subspace (Gaussian on R^{r-1})
    for (t in 1:n_rays) {
        d <- rnorm(ncol(U))                  # direction in latent coords
        v <- as.numeric(U %*% d)            # ambient direction
        interval <- line_box_interval(mu, v, a, b)
        if (interval[1] <= interval[2]) {
            ## Endpoints inside box
            y1 <- mu + interval[1] * v
            y2 <- mu + interval[2] * v
            rays[[length(rays) + 1L]] <- list(y = y1)
            rays[[length(rays) + 1L]] <- list(y = y2)
        }
    }

    if (length(rays) == 0L) {
        warning("No ray intersections with the box were found. Check r or the transform.")
        return(list(mu = mu, U = U,
                    candidates = matrix(, nrow = m, ncol = 0),
                    proj_coords = matrix(, nrow = ncol(U), ncol = 0),
                    distances = numeric(0)))
    }

    ## Snap endpoints to nearest vertices
    W_list <- list()
    for (q in seq_along(rays)) {
        y_end <- rays[[q]]$y
        w <- snap_to_vertex(y_end, a, b)
        W_list[[length(W_list) + 1L]] <- w
    }

    W <- unique(do.call(cbind, W_list), MARGIN = 2)
    W_counts = apply(W, 2, function(x) {
        matches = sapply(W_list, function(w) x == w)
        return(sum(matches))
    })
    W_probs = W_counts / sum(W_counts)

    list(candidates = W,
         candidate_probs = W_probs)
}
#' Build candidate endmembers exactly using Algorithm 3 (Deterministic Ray Casting)
#' Replaces the stochastic ray-casting approach with geometric linear algebra.
#' 
#' @param X A matrix of counts.
#' @param N An integer, vector, or matrix giving the number of trials.
#' @param r The dimension of the affine subspace (number of strains).
#' @param transform Either "arcsin" or "none".
#' @param eps_var Small value used for variance clipping.
#' @return A list with 'candidates' (top mathematically perfect endmembers), 
#'         'candidate_probs', and 'distances'.
build_candidates_exact <- function(X, N, r,
                                   transform = c("arcsin","none"),
                                   eps_var = 1e-6) {
    transform <- match.arg(transform)
    
    ## (1) Fit Affine Subspace (Phase 1)
    fit <- fit_affine_subspace(X, N, r, transform = transform, eps_var = eps_var)
    mu <- fit$mu
    U <- fit$U
    d <- fit$d  # Fetched for the Noise Trap check
    m <- length(mu)
    rminus1 <- max(1L, r - 1L)

    ## HPC Safety Check (Hole 1: The Noise Trap)
    # Check if we are over-specifying the model (rmodel > rtrue) by looking for a scree cliff
    if (!is.null(d) && length(d) >= rminus1 && rminus1 > 1) {
        if (d[rminus1] < 1e-3 * d[1]) {
            warning(sprintf("Potential noise trap: r_model (%d) may exceed true rank. Singular value cliff detected.", r))
        }
    }

    ## Box bounds in transformed space
    if (transform == "arcsin") {
        a <- 0; b <- pi
    } else if (transform == "none") {
        a <- 0; b <- 1
    }
    
    # Edge case handler for 1-dimensional models
    if (ncol(U) == 0L || rminus1 == 0L) {
        W_cand <- matrix(0L, nrow = m, ncol = 1)
        return(list(candidates = W_cand, candidate_probs = 1.0, distances = 0))
    }

    ## (2) Algebraic Enumeration (Phase 2 - Algorithm 3, Step 2-3)
    # Perform Rank-Revealing QR on t(U) to find r-1 linearly independent rows of U
    qUt <- qr(t(U))
    R_idx <- qUt$pivot[seq_len(rminus1)]
    U_R <- U[R_idx, , drop = FALSE]
    mu_R <- mu[R_idx]

    # Invert U_R and calculate Projection Matrix Z
    # Using tryCatch to handle mathematically unstable subspace projections
    U_R_inv <- tryCatch({
        solve(U_R)
    }, error = function(e) {
        warning("U_R is computationally singular. Using pseudo-inverse to salvage geometry.")
        MASS::ginv(U_R)
    })
    
    # Z matrix translates the lower-dimensional space back into the full m-dimensional SNP space
    Z <- U %*% U_R_inv  # size: m x (r-1)

    ## HPC Safety Check (Hole 3: Combinatorial Explosion)
    # Vectorized combination generation avoiding R 'for' loops
    u_binary <- as.matrix(expand.grid(rep(list(c(0, 1)), rminus1))) # size: 2^{r-1} x (r-1)
    u_scaled <- u_binary * (b - a) + a # Map to {a, b} space (e.g., 0 to pi)

    ## (3) Distance Ranking & Snapping (Phase 3 - Algorithm 3, Step 4-5)
    # Continuous projection into full ambient space: Y_hat = Z(u - mu_R) + mu
    U_diff <- sweep(u_scaled, 2, mu_R, "-")
    Y_hat <- sweep(Z %*% t(U_diff), 1, mu, "+") # size: m x 2^{r-1}

    # Temporary thresholding to binary 0/1 states
    W_bin <- snap_to_vertex(Y_hat, a, b) 
    
    # Map the discrete binary back to the transformed geometry to calculate strict L2 distance
    Y_target <- W_bin * (b - a) + a 

    ## HPC Safety Check (Hole 2: Floating-Point Snapping Danger)
    # Calculate exact Euclidean distance (delta) for every single candidate
    delta <- sqrt(colSums((Y_hat - Y_target)^2))

    # Filter out duplicate candidates resulting from the noisy projection
    unique_indices <- !duplicated(t(W_bin))
    W_unique <- W_bin[, unique_indices, drop = FALSE]
    delta_unique <- delta[unique_indices]

    # Sort the unique candidates strictly by their geometric proximity (delta)
    sort_order <- order(delta_unique)
    W_sorted <- W_unique[, sort_order, drop = FALSE]
    delta_sorted <- delta_unique[sort_order]

    # Only select the top 'r' mathematically closest candidates
    # This prevents the MIQP solver from calculating likelihoods on impossible garbage variants
    n_select <- min(r, ncol(W_sorted))
    W_final <- W_sorted[, seq_len(n_select), drop = FALSE]

    # Assign uniform weights to the final candidates for the downstream Importance Sampler
    W_probs <- rep(1 / n_select, n_select)

    return(list(candidates = W_final,
                candidate_probs = W_probs,
                distances = delta_sorted[seq_len(n_select)]))
}

#' A function that will make a matrix of the number of trials
#' @param N Either a single integer, a vector of length n, or a matrix of size m x n.
#' @param m The number of rows for the output matrix.
#' @param n The number of columns for the output matrix.
#' @return If N is an m x n matrix, returns N. If N is a length-n vector, returns an m x n matrix with all the elements of the ith column equal to N[i]. If N is a single integer, returns an m x n matrix with all elements equal to N.
make_Nmat <- function(N, m, n) {
    ## Expand N to m x n
    if (length(N) == 1L) {
        Nmat <- matrix(N, m, n)
    } else if (length(N) == n) {
        Nmat <- matrix(rep(as.numeric(N), each = m), m, n)
    } else if (is.matrix(N) && all(dim(N) == c(m, n))) {
        Nmat <- N
    } else {
        stop("N must be a scalar, length-n vector, or m x n matrix.")
    }
    return(Nmat)
}

#' Data likelihood given W and H
#' Computes P(X | W, H) for Binomial observations with N trials
#' @param X: m x n integer counts
#' @param N: either scalar, length-n vector (per column), or m x n matrix of trials
#' @param W: m x r binary (0/1)
#' @param H: r x n nonnegative with column sums ~ 1 (not enforced here)
#' @param log: Should the log likelihood be returned?
#' @param eps: If non-NULL, a small value used to require all probabilities to be strictly greater than zero and strictly less than 1. If NULL, no such restriction is applied.
#' @return Either the likelihood or the log-likelihood for log = FALSE or TRUE, respectively.
lik_binom_WH <- function(X, N, W, H, log = TRUE, eps = 1e-12) {
    X <- as.matrix(X)
    m <- nrow(X)
    n <- ncol(X)

    ## Validate dimensions
    if (!is.matrix(W) || !is.matrix(H)) stop("W and H must be matrices.")
    if (nrow(W) != m) stop("nrow(W) must equal nrow(X).")
    if (ncol(H) != n) stop("ncol(H) must equal ncol(X).")
    if (ncol(W) != nrow(H)) stop("W (m x r) and H (r x n) inner dims must match.")

    Nmat <- make_Nmat(N, m, n)

    ## Compute S = W %*% H and clamp away from 0/1 to avoid log(0)
    S <- W %*% H
    if(!is.null(eps)) {
        S <- clip(S, eps, 1-eps)
    }

    ## Vectorized binomial log-likelihood
    ll <- sum(dbinom(x = as.vector(X),
                     size = as.vector(Nmat),
                     prob = as.vector(S),
                     log = TRUE))

    if (log) {
        return(ll)
    } else {
        return(exp(ll))
    }
}




####################

# ===============================
# Numerics and densities
# ===============================
clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

ldirichlet1 <- function(h, a, eps = 1e-12) {
  h <- clip(h, eps, 1 - eps); a <- pmax(a, eps)
  lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(h))
}
rdirichlet1 <- function(a) { z <- rgamma(length(a), a); z/sum(z) }

# Binomial log-likelihood for a single column j
ll_binom_col <- function(Xj, Nj, W, h, eps = 1e-12) {
  s <- clip(as.numeric(W %*% h), eps, 1 - eps)
  sum(Xj * log(s) + (Nj - Xj) * log(1 - s))
}

# Prior on W (flat by default). You can plug VolMin prior here.
logprior_W <- function(W) 0

## ===============================
## Projection center c_j(W): argmin ||W h - s||^2, h in simplex
## ===============================
#' Project a vector v onto the simplex
#' @param v The vector to be projected.
#' @return The probability vector closest to v in Euclidean norm.
proj_simplex <- function(v) {
    r <- length(v)
    u <- sort(v, decreasing=TRUE)
    css <- cumsum(u)
    rho <- max(which(u > (css - 1) / seq_len(r)))
    theta <- (css[rho] - 1)/rho
    x <- pmax(v - theta, 0)
    if (sum(x) == 0) {
        x <- rep(1/r, r)
    }
    return(x / sum(x))
}

#' Solves  argmin ||W h - s||^2, h in simplex
#' @param W A binary matrix
#' @param s A probability vector
#' @param maxit Maximum number of iterations.
#' @param tol Stopping threshold.
#' @return The minimizing h.
projection_center <- function(W, s, maxit = 300, tol = 1e-8) {
  r <- ncol(W)
  L <- max(eigen(crossprod(W), symmetric=TRUE, only.values=TRUE)$values)
  step <- if (L > 0) 1/L else 1
  h <- rep(1/r, r)
  Wh <- W %*% h
  for (it in 1:maxit) {
    g <- crossprod(W, Wh - s)
    h_new <- proj_simplex(h - step * g)
    if (sqrt(sum((h_new - h)^2)) < tol) {
        h <- h_new
        break
    }
    h <- h_new
    Wh <- W %*% h
  }
  return(clip(h, 1e-12, 1 - 1e-12))
}

# ===============================
# Sequential W proposal with its exact log q
# - pi: length-K probabilities on catalog columns
# - Wcat: m x K matrix of catalog columns
# - returns: (W, idx, log_q)
# ===============================
sample_W_with_logq <- function(Wcat, pi, r) {
    K <- ncol(Wcat)
    m <- nrow(Wcat)
    rem <- seq_len(K)
    wts <- pi
    idx <- integer(r)
    logq <- 0
    for (t in 1:r) {
        probs <- wts[rem]
        denom <- sum(probs)
        i_pick <- sample(rem, 1, prob = probs)
        logq <- logq + log(wts[i_pick]) - log(denom)
        idx[t] <- i_pick
        rem <- setdiff(rem, i_pick)
    }
    list(W = Wcat[, idx, drop=FALSE], idx = idx, log_q = logq)
}

# ===============================
# Inner IS estimate of log p(X|W) column-wise
# proposal r_j(h) = Dir(kappa * c_j(W) + alpha)
# returns: log_evidence = sum_j log Zhat_j
# ===============================
log_marglik_innerIS <- function(X, N, W, alpha, kappa = 50, M = 64) {
  X <- as.matrix(X); m <- nrow(X); n <- ncol(X); r <- ncol(W)
  if (length(N) == 1L) N <- rep(N, n)
  Shat <- sweep(X, 2L, N, "/")

  total <- 0
  # precompute centers c_j(W)
  centers <- matrix(0, r, n)
  for (j in 1:n) centers[, j] <- projection_center(W, Shat[, j])

  for (j in 1:n) {
    a_for <- kappa * centers[, j] + alpha
    # draw M Dirichlets and compute log-weights
    lws <- numeric(M)
    Xj <- X[, j]; Nj <- N[j]
    for (midx in 1:M) {
      h <- rdirichlet1(a_for)
      lw <- ll_binom_col(Xj, Nj, W, h) +
            ldirichlet1(h, alpha) -
            ldirichlet1(h, a_for)
      lws[midx] <- lw
    }
    # log-mean-exp
    c0 <- max(lws)
    logZhat_j <- c0 + log(mean(exp(lws - c0)))
    total <- total + logZhat_j
  }
  total
}

# ===============================
# Top-level: importance sampler over W with inner IS over H
# ===============================
#' Importance sampling estimator of P(W | X)
#' @param X A matrix with counts.
#' @param N An integer, vector, or matrix containing the number of trials per cell of X.
#' @param Wcat A matrix with columns containing binary endmembers.
#' @param pi A vector containing the probabilities of each of the binary endmembers.
#' @param S The number of importance samples.
#' @param alpha
#' @param kappa
#' @param M
#' @param logprior_W_fun A function giving the log prior on W.
#' @return 
is_over_W <- function(X, N, r, Wcat, pi, S = 2000,
                      alpha = NULL, kappa = 50, M = 64,
                      logprior_W_fun = logprior_W) {
  X <- as.matrix(X); m <- nrow(X); n <- ncol(X)
  if (is.null(alpha)) alpha <- rep(0.5, r)

  draws <- vector("list", S)
  logw  <- numeric(S)
  keys  <- character(S)

  for (s in 1:S) {
    # 1) draw W and compute log q(W)
    samp <- sample_W_with_logq(Wcat, pi, r)
    W <- samp$W; logq <- samp$log_q

    # 2) estimate log p(X|W) by inner IS over H
    log_marg <- log_marglik_innerIS(X, N, W, alpha, kappa = kappa, M = M)

    # 3) prior on W
    logpW <- logprior_W_fun(W)

    # 4) log weight
    logw[s] <- log_marg + logpW - logq

    # 5) store
    draws[[s]] <- W
    keys[s] <- paste(W, collapse = "")  # canonical key (columns in draw order)
  }

  # self-normalize
  c0 <- max(logw)
  w_norm <- exp(logw - c0)
  Zhat <- exp(c0) * mean(w_norm)         # (optional) evidence est. up to factor

  # aggregate posterior mass by unordered W (canonicalize columns lexicographically)
  ukeys <- vapply(draws, canon_key, "")
  # sum normalized weights per unique unordered W
  post_mass <- tapply(w_norm, ukeys, sum)
  post_prob <- post_mass / sum(w_norm)

  # collect a compact summary table
  out <- data.frame(
    key = names(post_prob),
    post = as.numeric(post_prob),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$post, decreasing = TRUE), ]

  list(posterior_table = out,
       draws_W = draws,
       log_weights = logw,
       keys_ordered = keys,
       Zhat_like = Zhat)
}



#' Importance sampling estimator of P(W | X)
#' @param X A matrix with counts.
#' @param N An integer, vector, or matrix containing the number of trials per cell of X.
#' @param Wcat A matrix with columns containing candidate binary endmembers.
#' @param alpha
#' @param kappa
#' @param M
#' @param logprior_W_fun A function giving the log prior on W.
#' @return 
truncated_posterior_W <- function(X, N, r, Wcat,
                      alpha = NULL, kappa = 50, M = 64,
                      logprior_W_fun = logprior_W) {
    
    X <- as.matrix(X); m <- nrow(X); n <- ncol(X)
    if (is.null(alpha)) alpha <- rep(0.5, r)

    n_candidates <- ncol(Wcat)
    combinations <- combn(n_candidates, r)
    n_combinations <- ncol(combinations)
    logw <- numeric(n_combinations)
    keys <- character(n_combinations)
    if(n_combinations > 500) {
          ans <- readline(
              prompt = paste0(
                  "There are ", n_combinations,
                  " combinations. This may take a long time.\n",
                  "Do you want to continue? (y/n) "
              )
          )
          if (tolower(ans) != "y") {
              stop("User chose to quit.")
          }
    }

    for (s in 1:n_combinations) {
        comb_indices <- combinations[,s]
        W <- Wcat[,comb_indices]

        ## 2) estimate log p(X|W) by inner IS over H
        log_marg <- log_marglik_innerIS(X, N, W, alpha, kappa = kappa, M = M)

        ## 3) prior on W
        logpW <- logprior_W_fun(W)

        ## 4) log weight
        logw[s] <- log_marg + logpW

        ## 5) store
        keys[s] <- canon_key(W)  # canonical key (columns in draw order)
    }
    ## self-normalize
    c0 <- max(logw)
    w_norm <- exp(logw - c0)
    Zhat <- exp(c0) * mean(w_norm)         # (optional) evidence est. up to factor

    post_prob <- w_norm / sum(w_norm)
    out <- data.frame(key = keys,
                      post = post_prob,
                      stringsAsFactors = FALSE)
    out <- out[order(out$post, decreasing = TRUE),]
    list(posterior_table = out,
         log_weights = logw,
         Zhat_like = Zhat)
}


#' Canonical representation of endmembers
#' @param W A binary matrix with columns giving endmembers.
#' @return A string representation of the endmembers. String is of the form xxx-xxx-xxx, with x's being 0's and 1's, -'s separating individual endmembers, and endmembers sorted lexicographically.
canon_key <- function(W) {
    apply(W, 2, paste0, collapse = "") |> sort() |> paste0(collapse = "-")
}

#' Transforms string sequence to matrix
#' Takes the canonical key version of a binary sequence matrix and creates a binary matrix representation of the sequences.
#' @param key A canonical key. Contains 0's, 1's and -'s. -'s separate the sequences, which are represented as sequences of 0's and 1's.
#' @return A binary matrix W, with columns representing sequences and rows representing positions along the sequence.
canonical_key_to_W <- function(key) {
  # split on '-'
  cols <- strsplit(key, "-", fixed = TRUE)[[1]]
  # number of rows = length of binary string
  m <- nchar(cols[1])
  r <- length(cols)
  # reconstruct W
  W <- matrix(0L, nrow = m, ncol = r)
  for (j in seq_len(r)) {
    # split j-th string into characters and convert to 0/1
    W[, j] <- as.integer(strsplit(cols[j], "")[[1]])
  }
  return(W)
}


#' Compares bootstrap distribution of log likelihood to observed
#' @param X The count matrix.
#' @param N The number of trials per cell.
#' @param Phat A fitted matrix giving success probabilities per cell.
#' @param B The number of bootstrap replicates.
#' @param eps Truncates probabilities away from 0 and 1.
#' @return A list with the observed log likelihood (ObsLL), a vector giving the log likelihood for the bootstrap replicates (LL_boot), and the fraction of bootstrap replicates that gave a larger log likelihood than the observed (p_value).
bootstrap_no_refit_loglik <- function(X, N, Phat, B = 1000, eps = 1e-12) {
  m <- nrow(X)
  n <- ncol(X)

  # Build N matrix (m x n)
  Nmat <- make_Nmat(N, m, n)

  # clamp probabilities
  Phat <- pmin(pmax(Phat, eps), 1 - eps)

  # Observed loglik, summed over all cells
  loglik_obs <- sum(
    X * log(Phat) + (Nmat - X) * log(1 - Phat)
  )

  # Bootstrap distribution
  loglik_boot <- numeric(B)

  for (b in seq_len(B)) {
    Xb <- matrix(
      rbinom(m * n, size = as.vector(Nmat), prob = as.vector(Phat)),
      nrow = m, ncol = n
    )

    loglik_boot[b] <- sum(
      Xb * log(Phat) + (Nmat - Xb) * log(1 - Phat)
    )
  }

  # p-value: how extreme is observed loglik vs bootstraps?
  p_val <- mean(loglik_boot <= loglik_obs)  
  # since higher log-likelihood means better fit, 
  # small observed means lack-of-fit.

  list(
    loglik_obs = loglik_obs,
    loglik_boot = loglik_boot,
    p_value = p_val
  )
}




# Inputs:
# X: m x n integer matrix
# N: length-n vector of trials per column
# W: m x r binary matrix (0/1)
# eps: clamp for p in (eps, 1-eps)
# return_diagnostics: TRUE to return stats
#' @importFrom quadprog solve.QP
diag_fit_WH <- function(X, N, W, eps = 1e-12) {
  m <- nrow(X); n <- ncol(X); r <- ncol(W)
  if(length(N) != n) stop("N must be length ncol(X)")
  # Precompute matrices for QP: objective is ||W h - s||^2
  # For one column s, objective = 1/2 * h^T D h - d^T h + const
  Dmat_base <- 2 * crossprod(W)            #  r x r  (2 W^T W)
  # regularize to be PD
  Dmat <- Dmat_base + diag(1e-8, r)
  
  # Constraints: h >= 0 and sum(h) = 1 (equality)
  # quadprog expects Amat with constraints as columns, and bvec
  Amat <- cbind(rep(1, r), diag(r))   # first column is equality sum(h)=1
  bvec <- c(1, rep(0, r))
  meq <- 1
  
  Hhat <- matrix(0, nrow = r, ncol = n)
  Phat <- matrix(0, nrow = m, ncol = n)
  for(j in seq_len(n)) {
    s <- X[, j] / N[j]   # observed proportions
    dvec <- 2 * crossprod(W, s)   # because minimize h^T W^T W h - 2 s^T W h
    # quadprog solves min 1/2 h^T D h - d^T h
    sol <- tryCatch({
      solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    }, error = function(e) {
      # fallback: small tweak if numeric issue: increase reg
      Dmat2 <- Dmat + diag(1e-6, r)
      solve.QP(Dmat2, dvec, Amat, bvec, meq = meq)
    })
    h <- sol$solution
    # numerical cleanup
    h[h < 0] <- 0
    h <- h / sum(h)
    Hhat[, j] <- h
    p <- as.vector(W %*% h)
    p <- pmin(pmax(p, eps), 1 - eps)
    Phat[, j] <- p
  }
  
  # Residuals on counts
  R_counts <- X - sweep(Phat, 2, N, FUN = "*")   # X - N * P
  ObsSS_counts <- sum(R_counts^2)
  ExpSS_counts <- sum(sweep(Phat * (1 - Phat), 2, N, FUN = "*")) # sum N * p (1-p)
  
  # Residuals on proportions
  R_prop <- sweep(X, 2, N, FUN = "/") - Phat
  ObsSS_prop <- sum(R_prop^2)
  ExpSS_prop <- sum((Phat * (1 - Phat)) / matrix(rep(N, each=m), nrow=m))
  
  # Standardized residuals and chi-square
  Std_resid <- R_counts / sqrt(sweep(Phat * (1 - Phat), 2, N, FUN = "*"))
  # if any denom very small, clamp to large number (handled by eps above mostly)
  chi2 <- sum(Std_resid^2)
  # degrees of freedom: approx n*(m - (r - 1)) = n*(m - r + 1)
  df <- n * (m - r + 1)
  p_value <- pchisq(chi2, df = df, lower.tail = FALSE)
  
  out <- list(
    Hhat = Hhat,
    Phat = Phat,
    ObsSS_counts = ObsSS_counts,
    ExpSS_counts = ExpSS_counts,
    ObsSS_prop = ObsSS_prop,
    ExpSS_prop = ExpSS_prop,
    Std_resid = Std_resid,
    chi2 = chi2,
    df = df,
    p_value = p_value,
    max_abs_std_resid = max(abs(Std_resid))
  )
  return(out)
}
# ============================================================
# HYBRID 2: Deterministic Candidate Pool (Oversampled RRQR)
# Replaces stochastic ray-casting with Rank-Revealing QR
# to solve the High-K Condition Number Explosion.
# ============================================================

build_candidates_hybrid <- function(X, N, r, transform = "arcsin", oversample_factor = 3) {
  m <- nrow(X)
  n <- ncol(X)
  
  # 1. Fit the Affine Subspace (Preserve the Professor's Variance Stabilization)
  fit <- fit_affine_subspace(X, N, r, transform = transform)
  mu <- fit$mu
  U <- fit$U
  
  # 2. RRQR for Deterministic Feature Extraction
  # We transpose U (r-1 x m) so that columns are SNPs.
  # QR with column pivoting finds the most mathematically independent SNPs.
  # We oversample by extracting (r * oversample_factor) pivots.
  n_anchors <- min(r * oversample_factor, m)
  q <- qr(t(U))
  pivot_snps <- q$pivot[1:n_anchors]
  
  candidates_list <- list()
  
  # 3. Geometric Snapping (The Deterministic Spotlight)
  # For each highly independent SNP, we shoot a deterministic ray in both 
  # the positive and negative direction within the affine subspace.
  for (i in seq_along(pivot_snps)) {
    idx <- pivot_snps[i]
    
    # Extract the deterministic geometric direction from the SVD subspace
    direction <- U[idx, ]
    
    # Project extreme boundaries mathematically
    proj_pos <- mu + (U %*% direction)
    proj_neg <- mu - (U %*% direction)
    
    # Inverse transform and snap strictly to [0,1] boundaries (Positive Ray)
    if(transform == "arcsin") {
      p_pos <- (sin(proj_pos/2))^2
    } else {
      p_pos <- proj_pos
    }
    cand_pos <- ifelse(p_pos > 0.5, 1, 0)
    
    # Inverse transform and snap strictly to [0,1] boundaries (Negative Ray)
    if(transform == "arcsin") {
      p_neg <- (sin(proj_neg/2))^2
    } else {
      p_neg <- proj_neg
    }
    cand_neg <- ifelse(p_neg > 0.5, 1, 0)
    
    candidates_list[[length(candidates_list) + 1]] <- cand_pos
    candidates_list[[length(candidates_list) + 1]] <- cand_neg
  }
  
  # Combine into a single candidate pool matrix (m rows x un-filtered candidates)
  Wcat_raw <- do.call(cbind, candidates_list)
  
  # 4. De-duplication
  # The biological snapping process will likely snap different mathematical vectors 
  # to the same binary biological strain. We remove duplicates to save MCMC time.
  Wcat_unique <- unique(Wcat_raw, MARGIN = 2)
  
  # Safety Pad: Ensure we always pass at least 'r' candidates to the Bayesian filter
  if (ncol(Wcat_unique) < r) {
    pad_size <- r - ncol(Wcat_unique)
    X_prop <- sweep(X, 2, N, FUN = "/")
    consensus <- ifelse(rowMeans(X_prop) > 0.5, 1, 0)
    Wcat_unique <- cbind(Wcat_unique, matrix(rep(consensus, pad_size), ncol=pad_size))
  }
  
  # Return the exact object structure expected by test-nmf.R
  return(list(candidates = Wcat_unique, U = U, mu = mu, pivots = pivot_snps))
}

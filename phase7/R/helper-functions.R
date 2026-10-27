


#' Cube vertices in affine hull of columns of W
#'
#' @param W integer/numeric matrix n x k with entries 0/1
#' @return list with X (n x |V|), Z (n x |V|), and size |V|
#'         where X = aff(W) ∩ {0,1}^n.
cube_vertices_in_affine_hull <- function(W, tol = 1e-9, max_r = 25,
                                        return_as = c("matrix", "data.frame")) {
  return_as <- match.arg(return_as)
  W <- as.matrix(W)
  n <- nrow(W); k <- ncol(W)
  if (k == 0L) stop("W must have at least one column.")
  if (!all(W %in% c(0, 1))) stop("W must be 0/1.")

  # Choose basepoint w0 (first column)
  w0 <- as.integer(W[, 1L])

  if (k == 1L) {
    # aff(W) is just the single point w0
    X <- matrix(w0, nrow = n, ncol = 1L)
    colnames(X) <- "v1"
    if (return_as == "data.frame") return(list(X = as.data.frame(t(X)), size = 1L))
    return(list(X = X, size = 1L))
  }

  # Differences D = W_{-1} - w0
  D <- sweep(W[, -1L, drop = FALSE], 1L, w0, FUN = "-")

  # Sign flip rows where w0_i = 1
  s <- ifelse(w0 == 0L, 1, -1)            # length n
  Dtilde <- D * s                         # rowwise multiply (equiv to diag(s) %*% D)

  # Enumerate Z = span(Dtilde) ∩ {0,1}^n using your existing function
  Z <- cube_vertices_in_span(Dtilde, tol = tol, max_r = max_r, return_as = "matrix")

  # Map back: X = w0 + S Z  (rowwise: x_i = w0_i + s_i * z_i)
  X <- sweep(Z * s, 1L, w0, FUN = "+")
  storage.mode(X) <- "integer"

  # Safety check: should be 0/1
  if (!all(X %in% c(0L, 1L))) stop("Internal error: produced non-binary X.")

  if (return_as == "data.frame") {
    return(list(X = as.data.frame(t(X)), size = ncol(X)))
  }
  list(X = X, Z = Z, size = ncol(X))
}

# Convenience: just the size
size_aff_cube_intersection <- function(W, tol = 1e-9, max_r = 25) {
  cube_vertices_in_affine_hull(W, tol = tol, max_r = max_r, return_as = "matrix")$size
}


#' Computes the cube vertices in the column span of a matrix
#' #' @param A A full-column-rank matrix.
#' @return The elements of {0,1}^n which are in the span of the columns of A (equivalently, in the affine hull of 0, a_1,...,a_p).
#' @export
cube_vertices_in_span <- function(A, tol = 1e-9, max_r = 25, return_as = c("matrix", "data.frame")) {
  return_as <- match.arg(return_as)

  A <- as.matrix(A)
  n <- nrow(A)
  if (n == 0L) stop("A must have at least one row.")
  if (ncol(A) == 0L) {
    # Span is {0}. Only vertex in {0,1}^n is the origin.
    V <- matrix(0L, nrow = n, ncol = 1L)
    colnames(V) <- "v1"
    if (return_as == "data.frame") return(as.data.frame(t(V)))
    return(V)
  }

  # 1) Find rank r and a column-basis matrix B (n x r)
  qA <- qr(A)
  r <- qA$rank
  if (r == 0L) {
    V <- matrix(0L, nrow = n, ncol = 1L)
    colnames(V) <- "v1"
    if (return_as == "data.frame") return(as.data.frame(t(V)))
    return(V)
  }
  if (r > max_r) {
    stop(sprintf("rank(A) = %d is > max_r = %d; enumeration over 2^r is too large.", r, max_r))
  }

  piv_cols <- qA$pivot[seq_len(r)]
  B <- A[, piv_cols, drop = FALSE]  # basis for col(A)

  # 2) Choose row index set I of size r so that B_I is invertible.
  #    Use QR on t(B) to pick r independent columns of t(B), i.e. independent rows of B.
  qBt <- qr(t(B))
  piv_rows <- qBt$pivot[seq_len(r)]
  I <- sort(piv_rows)
  J <- setdiff(seq_len(n), I)

  BI <- B[I, , drop = FALSE]  # r x r
  # Quick sanity check (numerical)
  if (qr(BI)$rank < r) stop("Failed to find an invertible B_I (numerical issue). Try different tol or a different basis.")

  BJ <- B[J, , drop = FALSE]  # (n-r) x r

  # 3) Compute M = B_J %*% solve(B_I)
  invBI <- solve(BI)
  M <- BJ %*% invBI  # (n-r) x r

  # Helper: check whether a numeric vector is (approximately) 0/1 entrywise
  is_near_01 <- function(v) {
    # accept if each coordinate is within tol of 0 or 1
    all(pmin(abs(v), abs(v - 1)) <= tol)
  }

  # 4) Enumerate u in {0,1}^r, form x = (u, M u) and keep those in {0,1}^n
  num_u <- 2^r
  # We'll store as list of integer vectors (then bind at the end)
  verts <- vector("list", length = 0L)

  # Loop over all bit patterns (r <= max_r <= 25 by default, so intToBits is fine)
  for (k in 0:(num_u - 1L)) {
    bits <- as.integer(intToBits(k))[seq_len(r)]
    u <- matrix(bits, nrow = r, ncol = 1L)

    v <- as.vector(M %*% u)  # length n-r
    if (!is_near_01(v)) next

    # Round to exact 0/1 after acceptance
    v01 <- as.integer(abs(v - 1) < abs(v))
    x <- integer(n)
    x[I] <- bits
    x[J] <- v01

    # Final exact check (protect against rare tolerance misclassifications)
    if (all(x %in% c(0L, 1L))) {
      verts[[length(verts) + 1L]] <- x
    }
  }

  if (length(verts) == 0L) {
    # Should not happen (0 is always in the span), but keep it safe.
    V <- matrix(0L, nrow = n, ncol = 1L)
    colnames(V) <- "v1"
    if (return_as == "data.frame") return(as.data.frame(t(V)))
    return(V)
  }

  V <- do.call(cbind, lapply(verts, as.integer))  # n x (#verts)
  colnames(V) <- paste0("v", seq_len(ncol(V)))

  if (return_as == "data.frame") {
    # return one vertex per row (often nicer)
    return(as.data.frame(t(V)))
  }
  V
}


 in_convex_hull <- function(W, x, tol = 1e-9, verbose = FALSE) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required. Install it with install.packages('lpSolve')")
  }

  W <- as.matrix(W)
  x <- as.numeric(x)
  n <- nrow(W)
  m <- ncol(W)

  if (length(x) != n)
    stop("x must have length nrow(W).")
  if (m == 0L)
    stop("W must have at least one column.")

  ## Constraints:
  ## 1) W %*% lambda = x          (n equalities)
  ## 2) sum(lambda) = 1           (1 equality)
  ## 3) lambda >= 0               (m inequalities)

  Aeq <- rbind(W, rep(1, m))     # (n+1) x m
  beq <- c(x, 1)

  Aineq <- diag(m)               # m x m
  bineq <- rep(0, m)

  const.mat <- rbind(Aeq, Aineq)
  const.dir <- c(rep("=", n + 1),
                 rep(">=", m))
  const.rhs <- c(beq, bineq)

  res <- lpSolve::lp(
    direction    = "min",
    objective.in = rep(0, m),     # feasibility LP
    const.mat    = const.mat,
    const.dir    = const.dir,
    const.rhs    = const.rhs
  )

  if (verbose) {
    message("lp status: ", res$status,
            " (0 = feasible, 2 = infeasible)")
  }

  if (res$status != 0) {
    return(list(in_hull = FALSE, lambda = NULL))
  }

  lambda <- res$solution

  ## Numerical verification (important near boundaries)
  ok_nonneg <- min(lambda) >= -10 * tol
  ok_sum    <- abs(sum(lambda) - 1) <= 10 * tol
  ok_eq     <- max(abs(W %*% lambda - x)) <= 10 * tol

  in_hull <- ok_nonneg && ok_sum && ok_eq
  if (!in_hull) lambda <- NULL

  list(in_hull = in_hull, lambda = lambda)
}

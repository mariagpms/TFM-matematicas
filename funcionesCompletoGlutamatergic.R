#Funciones auxiliares
normalize_weights_to_n <- function(w) {
  w <- as.numeric(w)
  n <- length(w)
  if (n == 0) return(w)
  sum_w <- sum(w, na.rm = TRUE) 
  if (sum_w == 0) return(rep(1, n)) 
  w * (n / sum_w)
}
angdiff <- function(x) {
  (x + pi) %% (2 * pi) - pi
}

# ============================================================
# 1) Split A/B circular (pares/impares en orden de t)
# ============================================================
make_split_k <- function(t, k = 3, screen_mod = 1) {
  # k: tamaño del ciclo (3 => 1/3 screen, 2/3 test)
  # screen_mod: qué posición del ciclo va a screening (1..k)
  stopifnot(k >= 2, screen_mod >= 1, screen_mod <= k)
  
  ord <- order(as.numeric(t) %% (2*pi))
  pos <- ((seq_along(ord) - 1) %% k) + 1  # 1..k repetido
  
  A <- ord[pos == screen_mod]      # screening
  B <- ord[pos != screen_mod]      # testing (resto)
  
  list(A = A, B = B, ord = ord, k = k, screen_mod = screen_mod)
}

# ============================================================
# Circular von Mises kernel smoother (Nadaraya–Watson)
# ============================================================
smooth_circular_vm <- function(t, x, kappa = 8, w = NULL, eps = 1e-12) {
  t <- as.numeric(t) %% (2*pi)
  x <- as.numeric(x)
  n <- length(t)
  if (length(x) != n) stop("x and t must have same length.")
  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w)
  if (length(w) != n) stop("w and t must have same length.")
  
  Delta <- outer(t, t, function(a, b) angdiff(a - b))
  K <- exp(kappa * cos(Delta))  # unnormalized von Mises kernel
  
  # weights for smoothing
  KW <- K * matrix(w, nrow = n, ncol = n, byrow = TRUE)
  denom <- rowSums(KW) + eps
  
  # Nadaraya–Watson estimate at each observed time point
  muhat <- rowSums(KW * matrix(x, nrow = n, ncol = n, byrow = TRUE)) / denom
  muhat
}

# ============================================================
# Metrics for ONE gene
# ============================================================
smooth_metrics_one_gene <- function(t, x, w = NULL, kappa = 8, eps = 1e-12) {
  t <- as.numeric(t) %% (2*pi)
  x <- as.numeric(x)
  if (length(x) != length(t)) stop("x must have same length as t.")
  n <- length(t)
  if (is.null(w)) w <- rep(1, n)
  w <- normalize_weights_to_n(w)
  
  mu0 <- sum(w * x) / sum(w)
  wrss0 <- sum(w * (x - mu0)^2)
  
  muhat <- smooth_circular_vm(t, x, kappa = kappa, w = w, eps = eps)
  wrss_smooth <- sum(w * (x - muhat)^2)
  
  r2 <- 1 - wrss_smooth / (wrss0 + eps)
  r2 <- max(0, min(1, r2))
  
  signal <- max(0, wrss0 - wrss_smooth)
  noise  <- wrss_smooth
  snr <- signal / (noise + eps)
  
  list(
    wrss0 = wrss0,
    wrss_smooth = wrss_smooth,
    r2_smooth = r2,
    snr_smooth = snr
  )
}

# ============================================================
# Metrics for ALL genes in a matrix genes x samples
# ============================================================
smooth_metrics_all_genes <- function(X, t, w = NULL, kappa = 8, eps = 1e-12, verbose = TRUE) {
  t <- as.numeric(t) %% (2*pi)
  n <- length(t)
  if (is.null(w)) w <- rep(1, n)
  w <- normalize_weights_to_n(w)
  # X: genes x samples
  stopifnot(ncol(X) == n, n == length(w))
  
  genes <- rownames(X)
  if (is.null(genes)) genes <- paste0("g", seq_len(nrow(X)))
  
  out <- data.frame(
    gene = genes,
    wrss0 = NA_real_,
    wrss_smooth = NA_real_,
    r2_smooth = NA_real_,
    snr_smooth = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(X))) {
    if (verbose && (i %% 500 == 0)) cat("gene", i, "/", nrow(X), "\n")
    m <- smooth_metrics_one_gene(t, X[i, ], w, kappa = kappa, eps = eps)
    out$wrss0[i] <- m$wrss0
    out$wrss_smooth[i] <- m$wrss_smooth
    out$r2_smooth[i] <- m$r2_smooth
    out$snr_smooth[i] <- m$snr_smooth
  }
  
  out
}

# ==============================================
### Para elegir el corte
snr_knee_filter <- function(metrics, snr_col = "snr_smooth") {
  snr <- metrics[[snr_col]]
  snr <- snr[is.finite(snr)]
  if (length(snr) < 5) stop("Not enough finite SNR values to compute a knee.")
  
  snr_sorted <- sort(snr)
  x <- seq_along(snr_sorted)
  x <- (x - min(x)) / (max(x) - min(x))
  y <- (snr_sorted - min(snr_sorted)) / (max(snr_sorted) - min(snr_sorted) + 1e-12)
  
  dist <- abs((y[length(y)] - y[1]) * x -
                (x[length(x)] - x[1]) * y +
                x[length(x)] * y[1] -
                y[length(y)] * x[1]) /
    sqrt((y[length(y)] - y[1])^2 + (x[length(x)] - x[1])^2 + 1e-12)
  
  knee_idx <- which.max(dist)
  snr_knee <- snr_sorted[knee_idx]
  
  keep_snr <- metrics[[snr_col]] >= snr_knee
  keep_snr[!is.finite(metrics[[snr_col]])] <- FALSE
  
  list(
    snr_knee = snr_knee,
    knee_idx = knee_idx,
    keep = keep_snr,
    snr_sorted = snr_sorted,
    dist = dist
  )
}

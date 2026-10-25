# Предобработка: фильтры детекции, бланков и RSD в QC; проверка дрейфа (leave-one-QC-out);
# нормализация PQN / MSTUS; импутация ½ минимума и log2; признаки «есть/нет».

rsd <- function(M) apply(M, 2, function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))

detection <- function(X, meta, groups) {
  qc <- meta$kind == "QC"
  out <- data.frame(det_QC = colMeans(!is.na(X[qc, , drop = FALSE])), row.names = colnames(X))
  for (g in groups) out[[paste0("det_", g)]] <- colMeans(!is.na(X[meta$study & meta$group == g, , drop = FALSE]))
  out
}

filter_features <- function(X, meta, cfg) {
  p <- cfg$preprocess; groups <- cfg$design$groups
  det <- detection(X, meta, groups)
  qc <- meta$kind == "QC"
  keep <- det$det_QC >= p$min_qc_detect &
          apply(det[paste0("det_", groups)], 1, max) >= p$min_group_detect
  log_msg(sprintf("Фильтр детекции (QC ≥ %.0f%%, ≥ %.0f%% хотя бы в одной группе): %d → %d",
                  100 * p$min_qc_detect, 100 * p$min_group_detect, ncol(X), sum(keep)))
  bl <- meta$kind == "blank"
  det$blank_ratio <- NA_real_
  if (any(bl)) {
    B <- X[bl, , drop = FALSE]; B[is.na(B)] <- 0
    mb <- col_medians(B); mb[mb == 0] <- NA
    ratio <- col_medians(X[qc, , drop = FALSE]) / mb
    ratio[is.na(ratio)] <- Inf
    det$blank_ratio <- ratio
    kb <- ratio >= p$blank_ratio
    log_msg(sprintf("Фильтр по бланкам (%d бланков, QC/бланк ≥ %g): удалено %d", sum(bl), p$blank_ratio, sum(keep & !kb)))
    keep <- keep & kb
  } else {
    log_msg("Бланков в таблице нет: фильтр по бланкам не применяется (выполнен на этапе MS-DIAL)")
  }
  det$RSD_QC_raw <- rsd(X[qc, , drop = FALSE])
  k2 <- keep & !is.na(det$RSD_QC_raw) & det$RSD_QC_raw <= p$max_qc_rsd
  log_msg(sprintf("Фильтр RSD_QC ≤ %.0f%%: %d → %d. Медианный RSD QC: %.1f%% до, %.1f%% после",
                  100 * p$max_qc_rsd, sum(keep), sum(k2),
                  100 * stats::median(det$RSD_QC_raw[keep], na.rm = TRUE),
                  100 * stats::median(det$RSD_QC_raw[k2], na.rm = TRUE)))
  det$pass_detection <- keep
  det$kept <- k2
  list(X = X[, k2, drop = FALSE], table = det)
}

# Локально-линейная регрессия с трикубическими весами (как LOWESS, frac = 1, без итераций устойчивости)
lowess_at <- function(x, y, x0, frac = 1) {
  n <- length(x); k <- max(2, floor(frac * n + 1e-10))
  vapply(x0, function(t) {
    d <- abs(x - t); h <- sort(d)[k]
    w <- if (h > 0) pmax(0, 1 - (d / h)^3)^3 else as.numeric(d == 0)
    if (sum(w > 0) < 2) return(stats::weighted.mean(y, w + 1e-12))
    xm <- sum(w * x) / sum(w); ym <- sum(w * y) / sum(w)
    sxx <- sum(w * (x - xm)^2)
    if (sxx <= 0) return(ym)
    ym + sum(w * (x - xm) * (y - ym)) / sxx * (t - xm)
  }, 0)
}

trend <- function(xq, yq, xall) {
  fit <- lowess_at(xq, yq, xall)
  lo <- min(xq); hi <- max(xq)
  fit[xall < lo] <- lowess_at(xq, yq, lo); fit[xall > hi] <- lowess_at(xq, yq, hi)
  fit
}

# Leave-one-QC-out: ошибка предсказания внутренних QC без коррекции (медиана остальных QC) и с LOWESS
drift_loo <- function(X, order, qc_ids) {
  L <- log2(X[qc_ids, , drop = FALSE]); xq <- order[qc_ids]
  inner <- order(xq)[-c(1, length(xq))]
  e0 <- e1 <- numeric(0)
  for (j in seq_len(ncol(L))) {
    y <- L[, j]
    for (k in inner) {
      m <- !is.na(y); m[k] <- FALSE
      if (is.na(y[k]) || sum(m) < 4) next
      e0 <- c(e0, abs(y[k] - stats::median(y[m])))
      e1 <- c(e1, abs(y[k] - trend(xq[m], y[m], xq[k])))
    }
  }
  c(none = stats::median(e0), lowess = stats::median(e1))
}

drift_correct <- function(X, order, qc_ids) {
  L <- log2(X); xq <- order[qc_ids]; xall <- order[rownames(X)]
  for (j in seq_len(ncol(L))) {
    yq <- L[qc_ids, j]; m <- !is.na(yq)
    if (sum(m) >= 4) L[, j] <- L[, j] - trend(xq[m], yq[m], xall) + stats::median(yq[m])
  }
  2^L
}

# PQN: референс — медианный спектр QC; фактор пробы — медиана отношений к референсу
pqn_factors <- function(X, ref_ids) {
  ref <- col_medians(X[ref_ids, , drop = FALSE])
  row_medians(sweep(X, 2, ref, "/"))
}

# MSTUS: сумма признаков, обнаруженных во всех пробах (исследуемых и QC), относительно медианы QC
mstus_factors <- function(X, qc_ids) {
  full <- colSums(is.na(X)) == 0
  tot <- rowSums(X[, full, drop = FALSE])
  list(f = tot / stats::median(tot[qc_ids]), n = sum(full))
}

normalize <- function(X, meta, method) {
  qc <- rownames(X)[meta[rownames(X), "kind"] == "QC"]
  if (method == "PQN") { f <- pqn_factors(X, qc); txt <- "PQN (референс — медианный спектр QC)" }
  else if (method == "MSTUS") { m <- mstus_factors(X, qc); f <- m$f; txt <- sprintf("MSTUS (%d признаков во всех пробах)", m$n) }
  else { f <- stats::setNames(rep(1, nrow(X)), rownames(X)); txt <- "без нормализации" }
  list(X = X / f, f = f, info = txt)
}

# Импутация ½ минимума признака (по исследуемым пробам и QC), затем log2
impute_log <- function(Xn) {
  mn <- apply(Xn, 2, min, na.rm = TRUE) / 2
  idx <- which(is.na(Xn), arr.ind = TRUE)
  Xn[idx] <- mn[idx[, 2]]
  log2(Xn)
}

presence_absence <- function(det, groups, lo = 0.5, hi = 0.8) {
  D <- as.matrix(det[paste0("det_", groups)])
  apply(D, 1, min) < lo & apply(D, 1, max) >= hi
}

# Модуль 5 — заполнение пропусков; модуль 6 — коррекция сигнала (EigenMS) и логарифмирование.

# Заполнение пропусков случайным лесом по алгоритму missForest (Stekhoven, Bühlmann, 2012):
# начальное заполнение средним; на каждой итерации для каждого признака с пропусками лес обучается на
# наблюдённых значениях (предикторы — остальные признаки) и предсказывает пропуски; остановка, когда
# относительное изменение заполненных значений перестаёт уменьшаться (возвращается предыдущая итерация).
# Все признаки итерации обновляются по матрице предыдущей итерации, поэтому результат не зависит от
# числа потоков. Расчёт ведётся в шкале log2 (распределения интенсивностей сильно скошены).
rf_impute <- function(X, maxiter = 10, ntree = 100, cores = 1, seed = 1) {
  L <- log2(X); miss <- is.na(L); p <- ncol(L)
  vars <- which(colSums(miss) > 0)
  if (!length(vars)) return(list(X = X, iterations = 0, nrmse_oob = NA, diff = numeric(0)))
  vars <- vars[order(colSums(miss)[vars])]
  Li <- L; for (j in vars) Li[miss[, j], j] <- mean(L[, j], na.rm = TRUE)
  set.seed(seed); seeds <- matrix(sample.int(1e8, maxiter * p), maxiter)
  diff_prev <- Inf; diffs <- numeric(0); oob <- NA
  for (it in seq_len(maxiter)) {
    cur <- Li
    res <- par_lapply(vars, function(j) {
      set.seed(seeds[it, j]); obs <- !miss[, j]
      rf <- randomForest::randomForest(x = cur[obs, -j, drop = FALSE], y = cur[obs, j], ntree = ntree,
                                       mtry = max(1, floor(sqrt(p - 1))))
      list(pred = stats::predict(rf, cur[!obs, -j, drop = FALSE]), nmse = utils::tail(rf$mse, 1) / stats::var(cur[obs, j]))
    }, cores)
    new <- cur
    for (k in seq_along(vars)) new[miss[, vars[k]], vars[k]] <- res[[k]]$pred
    d <- sum((new[miss] - cur[miss])^2) / sum(new[miss]^2)
    diffs <- c(diffs, d)
    log_msg(sprintf("    итерация %d: относительное изменение %.2e", it, d))
    if (d >= diff_prev) break                     # изменение перестало уменьшаться: остаётся предыдущая итерация (Li)
    Li <- new; diff_prev <- d
    oob <- sqrt(mean(vapply(res, `[[`, 0, "nmse"), na.rm = TRUE))
  }
  list(X = 2^Li, iterations = length(diffs), nrmse_oob = oob, diff = diffs)
}

half_min_impute <- function(X) apply(X, 2, function(v) { v[is.na(v)] <- min(v, na.rm = TRUE) / 2; v })

module5_impute <- function(s) {
  out <- module_dir(s, 5); d <- load_module(s, 4); X <- d$X; m <- s$impute
  pm <- 100 * mean(is.na(X))
  log_msg(sprintf("  Пропусков: %.1f%% (%d значений); метод: %s", pm, sum(is.na(X)), m$method))
  info <- list(method = m$method, missing_pct = pm)
  if (m$method == "rf") {
    r <- rf_impute(X, m$maxiter, m$ntree, s$project$cores, s$project$seed + 5)
    Xi <- r$X; info$iterations <- r$iterations; info$nrmse_oob <- r$nrmse_oob
    log_msg(sprintf("  Случайный лес (алгоритм missForest): итераций %d, оценка NRMSE (OOB) %.3f", r$iterations, r$nrmse_oob))
  } else if (m$method == "missForest") {
    if (!requireNamespace("missForest", quietly = TRUE)) stop("impute$method = missForest, но пакет missForest не установлен (или выберите rf)")
    r <- missForest::missForest(as.data.frame(X), maxiter = m$maxiter, ntree = m$ntree)
    Xi <- as.matrix(r$ximp); dimnames(Xi) <- dimnames(X); info$nrmse_oob <- unname(r$OOBerror[1])
    log_msg(sprintf("  missForest: NRMSE (OOB) %.3f", info$nrmse_oob))
  } else Xi <- half_min_impute(X)
  stopifnot(!anyNA(Xi))
  write_csv(data.frame(sample = rownames(Xi), group = d$meta$group, as.data.frame(Xi), check.names = FALSE), file.path(out, "imputed.csv"))
  write_csv(as.data.frame(info), file.path(out, "impute_info.csv"))
  save_module(s, 5, list(X = Xi, info = d$info, meta = d$meta, missing = is.na(X), impute = info))
}

# EigenMS (Karpievitch et al., Bioinformatics 2009): остатки модели «признак ~ группа» раскладываются
# сингулярным разложением; число значимых трендов систематической ошибки определяется перестановочным
# тестом (Buja, Eyuboglu 1992; как в sva); тренды вычитаются из данных. Векторы трендов ортогональны
# дизайну групп, поэтому средние групп не меняются. ВНИМАНИЕ: метод использует метки групп — оценки
# точности классификации на скорректированных данных завышены (см. модуль 7, вложенная проверка).
eigenms <- function(L, groups, n_perm = 200, alpha = 0.05, seed = 1) {
  Y <- t(L); n <- ncol(Y)
  Xd <- stats::model.matrix(~ factor(groups))
  Q <- diag(n) - Xd %*% solve(crossprod(Xd)) %*% t(Xd)
  R <- Y %*% Q
  sv <- svd(R)
  kmax <- min(n - ncol(Xd), length(sv$d), 50)
  ds <- (sv$d^2 / sum(sv$d^2))[seq_len(kmax)]
  set.seed(seed); cnt <- numeric(kmax)
  for (b in seq_len(n_perm)) {
    Rp <- t(apply(R, 1, sample)) %*% Q
    d0 <- svd(Rp, nu = 0, nv = 0)$d
    cnt <- cnt + ((d0^2 / sum(d0^2))[seq_len(kmax)] >= ds)
  }
  p <- (cnt + 1) / (n_perm + 1)
  p <- cummax(p)                                  # последовательная проверка: тренд k значим, только если значимы 1..k-1
  k <- sum(p <= alpha)
  Yn <- if (k > 0) Y - Y %*% sv$v[, seq_len(k), drop = FALSE] %*% t(sv$v[, seq_len(k), drop = FALSE]) else Y
  list(L = t(Yn), n_trends = k, p = p, var_explained = ds, trends = sv$v[, seq_len(max(k, 1)), drop = FALSE])
}

pqn_factors <- function(X) {
  ref <- apply(X, 2, stats::median)
  apply(sweep(X, 2, ref, "/"), 1, stats::median)
}

sample_boxplot <- function(L, meta, title) {
  ids <- rownames(L)
  d <- data.frame(sample = factor(rep(ids, ncol(L)), levels = ids), value = as.vector(L), group = rep(meta[ids, "group"], ncol(L)))
  ggplot2::ggplot(d, ggplot2::aes(sample, value, fill = group)) + ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = group_colors(d$group)) +
    ggplot2::labs(x = "Образцы (в порядке ввода)", y = "Интенсивность (log2)", fill = NULL, title = title) + theme_um() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

module6_correct <- function(s) {
  out <- module_dir(s, 6); d <- load_module(s, 5); X <- d$X; meta <- d$meta; cc <- s$correct
  if (cc$normalization == "PQN") {
    f <- pqn_factors(X); X <- X / f
    log_msg(sprintf("  PQN: факторы разведения %.2f–%.2f", min(f), max(f)))
  }
  L0 <- if (isTRUE(cc$log2)) log2(X) else X
  info <- list(normalization = cc$normalization, log2 = isTRUE(cc$log2), method = cc$method)
  if (cc$method == "eigenms") {
    e <- eigenms(L0, meta$group, cc$eigen_perm, cc$eigen_alpha, s$project$seed + 6)
    L <- e$L; info$n_trends <- e$n_trends
    log_msg(sprintf("  EigenMS: найдено и удалено трендов систематической ошибки: %d (перестановок %d)", e$n_trends, cc$eigen_perm))
    write_csv(data.frame(trend = seq_along(e$p), var_explained = e$var_explained, p = e$p), file.path(out, "eigenms_trends.csv"))
  } else L <- L0
  save_fig(s, patchwork::wrap_plots(sample_boxplot(L0, meta, "(a) до коррекции"), sample_boxplot(L, meta, "(b) после коррекции"), ncol = 1),
           out, "distribution_by_sample", 24, 14)
  write_csv(data.frame(sample = rownames(L), group = meta$group, as.data.frame(L), check.names = FALSE), file.path(out, "corrected_log2.csv"))
  write_csv(as.data.frame(info), file.path(out, "correct_info.csv"))
  save_module(s, 6, list(L = L, L_raw = L0, info = d$info, meta = meta, correct = info))
}

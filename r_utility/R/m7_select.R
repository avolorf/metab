# Модуль 7. Статистика и отбор признаков (гибридный алгоритм этапа 1, исправленный):
# ОСФ (одномерный статистический фильтр) → 4 модели caret → частота встречаемости (СПФ) → РОП (RFE).
# Поправка на множественные сравнения — по всем признакам, метод — параметр (по умолчанию Холм).
# Вложенная кросс-валидация повторяет весь отбор внутри обучающих блоков и даёт несмещённую оценку.

# ------------------------------------------------------------------ ОСФ
uvf <- function(L, groups, p_adjust = "holm", alpha = 0.05) {
  g <- factor(groups); k <- nlevels(g)
  st1 <- p_adjust == "stage1"
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  res <- lapply(seq_len(ncol(L)), function(j) {
    v <- L[, j]
    p_norm <- if (st1) safe(stats::shapiro.test(v)$p.value)                       # как на этапе 1: по всей выборке
              else min(vapply(levels(g), function(l) { x <- v[g == l]; if (length(x) < 3 || stats::sd(x) == 0) NA_real_ else safe(stats::shapiro.test(x)$p.value) }, 0), na.rm = TRUE)
    if (!is.finite(p_norm)) p_norm <- 0
    p_hom <- safe(stats::bartlett.test(v, g)$p.value); if (is.na(p_hom)) p_hom <- 0
    test <- if (p_norm < 0.05) "Уилкоксон" else if (p_hom < 0.05) "Уэлч" else "Стьюдент"
    pin <- if (st1) "BH" else "holm"                                                 # поправка между парами групп (k > 2)
    p <- if (k == 2) {
      a <- v[g == levels(g)[1]]; b <- v[g == levels(g)[2]]
      switch(test, "Уилкоксон" = safe(stats::wilcox.test(a, b, exact = FALSE)$p.value),
             "Уэлч" = safe(stats::t.test(a, b, var.equal = FALSE)$p.value), "Стьюдент" = safe(stats::t.test(a, b, var.equal = TRUE)$p.value))
    } else {
      pm <- switch(test, "Уилкоксон" = safe(stats::pairwise.wilcox.test(v, g, p.adjust.method = pin, exact = FALSE)$p.value),
                   "Уэлч" = safe(stats::pairwise.t.test(v, g, p.adjust.method = pin, pool.sd = FALSE)$p.value),
                   "Стьюдент" = safe(stats::pairwise.t.test(v, g, p.adjust.method = pin, pool.sd = TRUE)$p.value))
      min(pm, na.rm = TRUE)
    }
    c(p_norm = p_norm, p_homog = p_hom, p = p, test = match(test, c("Уилкоксон", "Уэлч", "Стьюдент")))
  })
  r <- as.data.frame(do.call(rbind, res)); r$id <- colnames(L)
  r$test <- c("Уилкоксон", "Уэлч", "Стьюдент")[r$test]
  r$p_adj <- if (p_adjust %in% c("stage1", "none")) r$p else stats::p.adjust(r$p, p_adjust)
  r$selected <- !is.na(r$p_adj) & r$p_adj <= alpha
  if (k == 2) r$log2FC <- colMeans(L[g == levels(g)[2], , drop = FALSE]) - colMeans(L[g == levels(g)[1], , drop = FALSE])
  r[, c("id", setdiff(names(r), "id"))]
}

# ------------------------------------------------------------------ модели caret и частота встречаемости
nb_funcs <- function() {                       # наивный байесовский классификатор (e1071) для caret::rfe
  f <- caret::rfFuncs
  f$fit <- function(x, y, first, last, ...) e1071::naiveBayes(x, y)
  f$pred <- function(object, x) {
    pr <- stats::predict(object, x, type = "raw")
    data.frame(pred = factor(colnames(pr)[max.col(pr, ties.method = "first")], levels = colnames(pr)), pr, check.names = FALSE)
  }
  f$rank <- function(object, x, y) {
    v <- caret::filterVarImp(x, y); v$Overall <- rowMeans(v)
    v <- v[order(v$Overall, decreasing = TRUE), , drop = FALSE]; v$var <- rownames(v); v
  }
  f
}

train_models <- function(df, s, seed) {
  suppressPackageStartupMessages(library(caret))       # knn3, plsda и др. вызываются caret без префикса
  sl <- s$select
  B <- sl$cv_folds * sl$cv_repeats
  fits <- list()
  for (m in sl$models) {
    ctrl <- caret::trainControl(method = "repeatedcv", number = sl$cv_folds, repeats = sl$cv_repeats,
                                seeds = caret_seeds(B, max(sl$tune_length, 25), seed + match(m, sl$models)))
    fits[[m]] <- tryCatch(suppressWarnings(caret::train(Label ~ ., data = df, method = m, metric = "Accuracy",
                                                        trControl = ctrl, tuneLength = sl$tune_length)),
                          error = function(e) { log_msg("    модель ", m, " не обучена: ", conditionMessage(e)); NULL })
  }
  fits[!vapply(fits, is.null, TRUE)]
}

top_features <- function(fit, n) {
  imp <- caret::varImp(fit, scale = FALSE)$importance
  rownames(imp) <- gsub("`", "", rownames(imp))
  imp$sum <- rowSums(imp)
  utils::head(rownames(imp)[order(-imp$sum)], n)
}

run_rfe <- function(df, s, seed) {
  sl <- s$select; x <- df[, -1, drop = FALSE]; y <- df$Label
  sizes <- seq_len(ncol(x)); B <- sl$rfe_folds * sl$rfe_repeats
  set.seed(seed)
  seeds <- c(lapply(seq_len(B), function(i) sample.int(1e6, length(sizes) + 1)), list(sample.int(1e6, 1)))
  ctrl <- caret::rfeControl(functions = if (sl$rfe_method == "nb") nb_funcs() else caret::rfFuncs, method = "repeatedcv",
                            number = sl$rfe_folds, repeats = sl$rfe_repeats, seeds = seeds, verbose = FALSE)
  suppressWarnings(caret::rfe(x, y, sizes = sizes, metric = "Accuracy", rfeControl = ctrl))
}

# Весь отбор: ОСФ → модели → частота → РОП. Возвращает признаки каждого шага.
select_features <- function(L, groups, s, seed, verbose = TRUE) {
  sl <- s$select; say <- function(...) if (verbose) log_msg(...)
  u <- uvf(L, groups, sl$p_adjust, sl$alpha)
  f_uvf <- u$id[u$selected]
  say(sprintf("  ОСФ (поправка: %s, α = %g): %d из %d признаков", sl$p_adjust, sl$alpha, length(f_uvf), ncol(L)))
  out <- list(uvf = u, f_uvf = f_uvf, f_spf = f_uvf, f_rfe = f_uvf, fits = list(), freq = NULL, rfe = NULL)
  if (length(f_uvf) < 2) { say("  Меньше двух признаков после ОСФ: модели и РОП пропущены"); return(out) }
  df <- data.frame(Label = factor(groups), L[, f_uvf, drop = FALSE], check.names = FALSE)
  out$fits <- train_models(df, s, seed)
  tops <- lapply(out$fits, top_features, n = sl$top_n)
  freq <- sort(table(unlist(tops)), decreasing = TRUE)                  # частота встречаемости (вместо tuple::tuplicated)
  out$freq <- data.frame(id = names(freq), n_models = as.integer(freq), stringsAsFactors = FALSE)
  out$f_spf <- out$freq$id[out$freq$n_models >= min(sl$min_freq, length(out$fits))]
  say(sprintf("  СПФ (в top-%d не менее чем у %d из %d моделей): %d признаков", sl$top_n, sl$min_freq, length(out$fits), length(out$f_spf)))
  if (length(out$f_spf) < 2) { out$f_rfe <- out$f_spf; return(out) }
  out$rfe <- run_rfe(data.frame(Label = factor(groups), L[, out$f_spf, drop = FALSE], check.names = FALSE), s, seed + 100)
  out$f_rfe <- caret::predictors(out$rfe)                              # ds_rfe собирается из результата RFE
  say(sprintf("  РОП (%s): %d признаков, точность CV %.3f", if (sl$rfe_method == "nb") "наивный Байес" else "случайный лес",
              length(out$f_rfe), max(out$rfe$results$Accuracy)))
  out
}

mean_auroc <- function(L, groups) {           # средняя AUROC отдельных признаков (описательная, как на этапе 1)
  v <- caret::filterVarImp(as.data.frame(L), factor(groups))
  mean(rowMeans(v))
}

auc_binary <- function(y, score) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Вложенная кросс-валидация: отбор признаков и классификатор обучаются только на обучающем блоке
nested_cv <- function(L, groups, s, seed) {
  sl <- s$select; g <- factor(groups); lev <- levels(g)
  pred <- list(); chosen <- list(); k <- 0
  for (r in seq_len(sl$nested_repeats)) {
    set.seed(seed + r)
    folds <- integer(length(g))
    for (l in lev) { i <- which(g == l); folds[i[sample.int(length(i))]] <- rep_len(seq_len(sl$nested_folds), length(i)) }
    for (f in seq_len(sl$nested_folds)) {
      k <- k + 1; te <- folds == f; tr <- !te
      sel <- select_features(L[tr, , drop = FALSE], g[tr], s, seed + 1000 * r + f, verbose = FALSE)
      feats <- sel$f_rfe; chosen[[k]] <- feats
      prob <- if (length(feats) == 0) matrix(rep(as.numeric(table(g[tr]) / sum(tr)), each = sum(te)), sum(te), dimnames = list(NULL, lev))
      else if (sl$rfe_method == "nb") stats::predict(e1071::naiveBayes(L[tr, feats, drop = FALSE], g[tr]), L[te, feats, drop = FALSE], type = "raw")
      else { set.seed(seed + k); stats::predict(randomForest::randomForest(L[tr, feats, drop = FALSE], g[tr]), L[te, feats, drop = FALSE], type = "prob") }
      pred[[k]] <- data.frame(repeat_ = r, fold = f, sample = rownames(L)[te], truth = g[te],
                              predicted = lev[max.col(prob, ties.method = "first")], prob, check.names = FALSE)
      log_msg(sprintf("    вложенная CV: повтор %d, блок %d — признаков %d", r, f, length(feats)))
    }
  }
  P <- do.call(rbind, pred)
  per_rep <- lapply(split(P, P$repeat_), function(d) c(accuracy = mean(d$predicted == d$truth),
                                                     auc = if (length(lev) == 2) auc_binary(as.integer(d$truth == lev[2]), d[[lev[2]]]) else NA))
  M <- do.call(rbind, per_rep)
  fr <- sort(table(unlist(chosen)), decreasing = TRUE)
  list(predictions = P,
       summary = data.frame(outer_folds = sl$nested_folds, repeats = sl$nested_repeats, input = sl$nested_input,
                            accuracy_mean = mean(M[, "accuracy"]), accuracy_sd = if (nrow(M) > 1) stats::sd(M[, "accuracy"]) else NA,
                            auc_mean = mean(M[, "auc"]), auc_sd = if (nrow(M) > 1) stats::sd(M[, "auc"]) else NA,
                            mean_n_features = mean(lengths(chosen))),
       frequency = data.frame(id = names(fr), share_of_folds = as.numeric(fr) / k, stringsAsFactors = FALSE))
}

# Необязательно: модель с ковариатами (пол, возраст), если они есть в метаданных и групп две
covariate_model <- function(L, meta) {
  g <- factor(meta$group); if (nlevels(g) != 2 || is.null(meta$sex) || is.null(meta$age)) return(NULL)
  D <- data.frame(group = as.numeric(g == levels(g)[2]), sex = as.numeric(meta$sex == "M"), age = meta$age)
  ok <- stats::complete.cases(D); if (sum(ok) < 10) return(NULL)
  D <- D[ok, ]; used <- names(D)[vapply(D, stats::sd, 0) > 0]
  Xd <- cbind(1, as.matrix(D[, used, drop = FALSE])); Y <- L[ok, , drop = FALSE]
  XtX <- solve(crossprod(Xd)); b <- XtX %*% crossprod(Xd, Y); res <- Y - Xd %*% b; dof <- nrow(Xd) - ncol(Xd)
  se <- sqrt(colSums(res^2) / dof * XtX[2, 2]); p <- 2 * stats::pt(abs(b[2, ] / se), dof, lower.tail = FALSE)
  data.frame(id = colnames(L), log2FC_adjusted = b[2, ], p = p, p_BH = stats::p.adjust(p, "BH"), covariates = paste(used[-1], collapse = "+"))
}

module7_select <- function(s) {
  out <- module_dir(s, 7); d <- load_module(s, 6); sl <- s$select
  meta <- d$meta; g <- factor(meta$group); info <- d$info
  ann <- function(df) merge(info[, c("id", "mz", "rt_min", "name")], df, by = "id", sort = FALSE)
  cl <- register_cores(s$project$cores); on.exit(stop_cores(cl))
  sel <- select_features(d$L, g, s, s$project$seed + 7)
  write_csv(ann(sel$uvf[order(sel$uvf$p), ]), file.path(out, "uvf_stats.csv"))
  tab <- function(f) data.frame(sample = rownames(d$L), Label = g, d$L[, f, drop = FALSE], check.names = FALSE)
  write_csv(tab(sel$f_uvf), file.path(out, "ds.csv")); write_csv(tab(sel$f_spf), file.path(out, "ds_d.csv"))
  write_csv(tab(sel$f_rfe), file.path(out, "ds_rfe.csv"))
  if (!is.null(sel$freq)) write_csv(ann(sel$freq), file.path(out, "feature_frequency.csv"))
  if (!is.null(sel$rfe)) write_csv(sel$rfe$results, file.path(out, "rfe_profile.csv"))
  acc <- data.frame()
  addacc <- function(fits, set) for (m in names(fits))
    acc <<- rbind(acc, data.frame(dataset = set, model = m, n_features = length(fits[[m]]$coefnames),
                                  accuracy_cv = mean(fits[[m]]$resample$Accuracy), accuracy_sd = stats::sd(fits[[m]]$resample$Accuracy)))
  addacc(sel$fits, "ds (после ОСФ)")
  if (isTRUE(sl$compare_stages)) {
    log_msg("  Точность моделей на всех признаках до и после коррекции (сравнение этапов)…")
    addacc(train_models(data.frame(Label = g, d$L_raw, check.names = FALSE), s, s$project$seed + 71), "ds_raw (до коррекции)")
    addacc(train_models(data.frame(Label = g, d$L, check.names = FALSE), s, s$project$seed + 72), "dsr (после коррекции)")
  }
  if (nrow(acc)) { write_csv(acc, file.path(out, "model_accuracy.csv"))
    for (i in seq_len(nrow(acc))) log_msg(sprintf("    %s, %s: точность CV %.3f", acc$dataset[i], acc$model[i], acc$accuracy_cv[i])) }
  steps <- data.frame(dataset = c("ds_raw (до коррекции)", "dsr (после коррекции)", "ds (после ОСФ)", "ds_d (после СПФ)", "ds_rfe (после РОП)"),
                      n_features = c(ncol(d$L_raw), ncol(d$L), length(sel$f_uvf), length(sel$f_spf), length(sel$f_rfe)))
  steps$mean_auroc_descriptive <- c(mean_auroc(d$L_raw, g), mean_auroc(d$L, g),
                                    vapply(list(sel$f_uvf, sel$f_spf, sel$f_rfe), function(f) if (length(f)) mean_auroc(d$L[, f, drop = FALSE], g) else NA, 0))
  write_csv(steps, file.path(out, "auroc_steps.csv"))
  log_msg("  Число признаков и средняя AUROC по этапам (описательная): ",
          paste(sprintf("%d (%.3f)", steps$n_features, steps$mean_auroc_descriptive), collapse = " → "))
  nest <- NULL
  if (isTRUE(sl$nested_cv)) {
    log_msg(sprintf("  Вложенная кросс-валидация (%d блоков × %d повтора, данные: %s)…", sl$nested_folds, sl$nested_repeats,
                    if (sl$nested_input == "imputed") "после заполнения пропусков, без EigenMS" else "после EigenMS — оценка завышена"))
    nest <- nested_cv(if (sl$nested_input == "imputed") d$L_raw else d$L, g, s, s$project$seed + 700)
    write_csv(nest$summary, file.path(out, "nested_cv.csv")); write_csv(nest$predictions, file.path(out, "nested_cv_predictions.csv"))
    write_csv(ann(nest$frequency), file.path(out, "nested_cv_feature_frequency.csv"))
    log_msg(sprintf("  Вложенная CV: точность %.3f, AUC %.3f (несмещённая оценка всей процедуры отбора)", nest$summary$accuracy_mean, nest$summary$auc_mean))
  }
  cov <- covariate_model(d$L_raw, meta)
  if (!is.null(cov)) { write_csv(ann(cov[order(cov$p), ]), file.path(out, "covariate_model.csv"))
    log_msg(sprintf("  Модель с поправкой на %s: q(BH) < 0,05 у %d признаков", cov$covariates[1], sum(cov$p_BH < 0.05))) }
  save_module(s, 7, list(L = d$L, L_raw = d$L_raw, info = info, meta = meta, sel = sel[c("uvf", "f_uvf", "f_spf", "f_rfe", "freq")],
                         rfe_results = if (!is.null(sel$rfe)) sel$rfe$results, accuracy = acc, steps = steps, nested = nest, covariate = cov))
}

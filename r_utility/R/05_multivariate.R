# Многомерный анализ без утечки информации: повторная стратифицированная кросс-валидация,
# отбор признаков внутри обучающего фолда, перестановочный тест, PCA.

FPR_GRID <- seq(0, 1, by = 0.01)

auc_score <- function(y, s) {
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_points <- function(y, s) {
  o <- order(-s); ys <- y[o]; ss <- s[o]
  last <- c(diff(ss) != 0, TRUE)                    # одна точка на каждый уникальный порог
  tp <- cumsum(ys == 1)[last]; fp <- cumsum(ys == 0)[last]
  list(fpr = c(0, fp / sum(y == 0)), tpr = c(0, tp / sum(y == 1)))
}

mean_roc <- function(y, s) { r <- roc_points(y, s); stats::approx(r$fpr, r$tpr, FPR_GRID, ties = "ordered", rule = 2)$y }

pareto_fit <- function(X) {
  mu <- colMeans(X); sdv <- sqrt(apply(X, 2, stats::sd)); sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(mu = mu, sd = sdv)
}

# Стратифицированное разбиение на k блоков
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) { i <- which(y == cl); i <- i[sample.int(length(i))]; f[i] <- rep_len(seq_len(k), length(i)) }
  f
}

# Логистическая регрессия с L2-штрафом (метод Ньютона; штраф 1/(2C)·||w||², свободный член без штрафа)
fit_logreg <- function(X, y, C = 1, iter = 100) {
  X1 <- cbind(1, X); w <- rep(0, ncol(X1)); P <- diag(c(0, rep(1 / C, ncol(X))), ncol(X1))
  for (it in seq_len(iter)) {
    p <- 1 / (1 + exp(-drop(X1 %*% w)))
    g <- crossprod(X1, p - y) + P %*% w
    H <- crossprod(X1 * (p * (1 - p)), X1) + P
    step <- solve(H, g); w <- w - drop(step)
    if (max(abs(step)) < 1e-8) break
  }
  w
}
predict_logreg <- function(w, X) 1 / (1 + exp(-drop(cbind(1, X) %*% w)))

residualize <- function(Xtr, Xte, ztr, zte) {
  zc <- ztr - mean(ztr)
  beta <- drop(crossprod(zc, sweep(Xtr, 2, colMeans(Xtr)))) / sum(zc^2)
  list(tr = Xtr - outer(zc, beta), te = Xte - outer(zte - mean(ztr), beta))
}

mw_p_fast <- function(Xtr, ytr) mw_matrix(Xtr[ytr == 1, , drop = FALSE], Xtr[ytr == 0, , drop = FALSE])$p

# Одна k-блочная CV. method: panel | panel_adj | plsda | covariate
cv_predict <- function(X, y, method, k, panel_k, z = NULL) {
  folds <- strat_folds(y, k); pred <- numeric(length(y)); chosen <- integer(0)
  for (f in seq_len(k)) {
    te <- folds == f; tr <- !te
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    if (method == "panel_adj") { r <- residualize(Xtr, Xte, z[tr], z[te]); Xtr <- r$tr; Xte <- r$te }
    if (method == "covariate") {
      mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, function(v) sqrt(mean((v - mean(v))^2))) + 1e-9
      w <- fit_logreg(sweep(sweep(Xtr, 2, mu), 2, sdv, "/"), y[tr])
      pred[te] <- predict_logreg(w, sweep(sweep(Xte, 2, mu), 2, sdv, "/"))
    } else if (method %in% c("panel", "panel_adj")) {
      sel <- order(mw_p_fast(Xtr, y[tr]))[seq_len(min(panel_k, ncol(X)))]
      chosen <- c(chosen, sel)
      A <- Xtr[, sel, drop = FALSE]
      mu <- colMeans(A); sdv <- apply(A, 2, function(v) sqrt(mean((v - mean(v))^2))) + 1e-9
      w <- fit_logreg(sweep(sweep(A, 2, mu), 2, sdv, "/"), y[tr])
      pred[te] <- predict_logreg(w, sweep(sweep(Xte[, sel, drop = FALSE], 2, mu), 2, sdv, "/"))
    } else {
      pf <- pareto_fit(Xtr)
      S <- sweep(sweep(Xtr, 2, pf$mu), 2, pf$sd, "/")
      d <- data.frame(y = y[tr]); d$X <- S
      m <- pls::plsr(y ~ X, ncomp = 2, data = d, scale = FALSE,
                     method = if (ncol(S) > nrow(S)) "widekernelpls" else "kernelpls")
      nd <- data.frame(row = seq_len(sum(te))); nd$X <- sweep(sweep(Xte, 2, pf$mu), 2, pf$sd, "/")
      pred[te] <- drop(stats::predict(m, newdata = nd, ncomp = 2))
    }
  }
  list(pred = pred, chosen = chosen)
}

repeated_cv <- function(X, y, method, seed, cfg, z = NULL) {
  mv <- cfg$multivariate
  set.seed(seed)
  aucs <- numeric(mv$cv_repeats); tprs <- matrix(0, mv$cv_repeats, length(FPR_GRID)); chosen <- integer(0); q2 <- 0
  for (r in seq_len(mv$cv_repeats)) {
    cv <- cv_predict(X, y, method, mv$cv_folds, mv$panel_k, z)
    aucs[r] <- auc_score(y, cv$pred); tprs[r, ] <- mean_roc(y, cv$pred); chosen <- c(chosen, cv$chosen)
    q2 <- q2 + (1 - sum((y - cv$pred)^2) / sum((y - mean(y))^2)) / mv$cv_repeats
  }
  list(auc = aucs, tpr = colMeans(tprs), chosen = chosen, q2 = q2)
}

# Перестановочный тест: метки перемешиваются, вся процедура CV (включая отбор признаков) повторяется
perm_test <- function(X, y, method, observed, seed, cfg, z = NULL) {
  n <- cfg$multivariate$n_perm
  if (!n) return(NA_real_)
  set.seed(seed); seeds <- sample.int(.Machine$integer.max, n)
  one <- function(s) {
    set.seed(s); yp <- sample(y)
    auc_score(yp, cv_predict(X, yp, method, cfg$multivariate$cv_folds, cfg$multivariate$panel_k, z)$pred)
  }
  null <- unlist(par_lapply(seeds, one, cores = cfg$cores))
  (sum(null >= observed) + 1) / (n + 1)
}

evaluate_models <- function(L, meta, case, ref, pool, label, seed, cfg, names = NULL) {
  idx <- c(case, ref); y <- c(rep(1, length(case)), rep(0, length(ref)))
  X <- L[idx, pool, drop = FALSE]
  z <- if (!is.null(meta$contam_index)) meta[idx, "contam_index"] else NULL
  methods <- c(panel = sprintf("панель top-%d + логистическая регрессия (отбор внутри CV)", cfg$multivariate$panel_k),
               panel_adj = sprintf("панель top-%d, поправка на индекс загрязнения внутри фолда", cfg$multivariate$panel_k),
               plsda = "PLS-DA, 2 компоненты (все соединения без экспозиций)")
  if (is.null(z)) methods <- methods[names(methods) != "panel_adj"]
  rows <- list(); curves <- list()
  for (m in names(methods)) {
    cv <- repeated_cv(X, y, m, seed, cfg, z)
    p <- perm_test(X, y, m, mean(cv$auc), seed + 1, cfg, z)
    sel <- ""
    if (length(cv$chosen)) {
      vc <- sort(table(cv$chosen), decreasing = TRUE)[seq_len(min(10, length(unique(cv$chosen))))]
      nm <- pool[as.integer(names(vc))]; if (!is.null(names)) nm <- names[nm]
      sel <- paste(sprintf("%s (%.0f%%)", nm, 100 * as.numeric(vc) / (cfg$multivariate$cv_repeats * cfg$multivariate$cv_folds)), collapse = "; ")
    }
    rows[[m]] <- data.frame(set = label, model = methods[[m]], n_case = length(case), n_ref = length(ref),
                            n_features = length(pool), AUC_mean = mean(cv$auc), AUC_sd = stats::sd(cv$auc),
                            Q2 = cv$q2, perm_p = p, n_perm = cfg$multivariate$n_perm, most_selected = sel)
    curves[[paste(label, m)]] <- list(tpr = cv$tpr, auc = mean(cv$auc), label = paste(label, methods[[m]]))
  }
  for (cv_name in intersect(c("sex", "contam_index"), names(meta))) {
    v <- meta[idx, cv_name]; v <- if (cv_name == "sex") c(M = 1, F = 0)[v] else as.numeric(v)
    ok <- !is.na(v)
    if (sum(ok) > 5 && stats::sd(v[ok]) > 0 && min(sum(y[ok]), sum(1 - y[ok])) >= cfg$multivariate$cv_folds) {
      Xv <- matrix(v[ok], ncol = 1)
      cv <- repeated_cv(Xv, y[ok], "covariate", seed, cfg)
      p <- perm_test(Xv, y[ok], "covariate", mean(cv$auc), seed + 1, cfg)
      rows[[cv_name]] <- data.frame(set = label, model = paste("только", cv_name, "(базовая модель, та же CV)"),
                                    n_case = sum(y[ok]), n_ref = sum(1 - y[ok]), n_features = 1,
                                    AUC_mean = mean(cv$auc), AUC_sd = stats::sd(cv$auc), Q2 = cv$q2, perm_p = p,
                                    n_perm = cfg$multivariate$n_perm, most_selected = "")
      curves[[paste(label, cv_name)]] <- list(tpr = cv$tpr, auc = mean(cv$auc), label = paste(label, "только", cv_name))
    }
  }
  list(table = do.call(rbind, rows), curves = curves)
}

pca_scores <- function(L, n = 3) {
  pf <- pareto_fit(L)
  p <- stats::prcomp(sweep(sweep(L, 2, pf$mu), 2, pf$sd, "/"), center = FALSE, scale. = FALSE, rank. = n)
  list(scores = p$x[, seq_len(n), drop = FALSE], var = (p$sdev^2 / sum(p$sdev^2))[seq_len(n)])
}

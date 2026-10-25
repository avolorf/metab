# Одномерная статистика: модель с ковариатами, подбор пар по возрасту, проверки устойчивости,
# правило «находки», демография.

bh <- function(p) stats::p.adjust(p, method = "BH")

comparison_sets <- function(meta, cm) {
  s <- meta[meta$study, ]
  if (!is.null(cm$sex)) s <- s[s$sex == cm$sex, ]
  list(case = s$sample[s$group == cm$case], ref = s$sample[s$group == cm$ref])
}

# Критерий Манна–Уитни (двусторонний, нормальное приближение с поправкой на связи и на непрерывность —
# как wilcox.test(exact = FALSE, correct = TRUE)) для всех столбцов сразу; AUC = U / (n1·n2).
mw_matrix <- function(A, B) {
  n1 <- nrow(A); n2 <- nrow(B); p <- ncol(A)
  if (n1 < 2 || n2 < 2) return(list(p = rep(NA_real_, p), auc = rep(NA_real_, p), fc = rep(NA_real_, p)))
  M <- rbind(A, B); n <- n1 + n2
  pv <- auc <- numeric(p)
  for (j in seq_len(p)) {
    r <- rank(M[, j]); U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
    t <- tabulate(match(M[, j], unique(M[, j])))
    s <- sqrt(n1 * n2 / 12 * ((n + 1) - sum(t^3 - t) / (n * (n - 1))))
    z <- (abs(U - n1 * n2 / 2) - 0.5) / s
    pv[j] <- if (s > 0) min(1, 2 * stats::pnorm(max(z, 0), lower.tail = FALSE)) else 1
    auc[j] <- U / (n1 * n2)
  }
  list(p = pv, auc = auc, fc = colMeans(A) - colMeans(B))
}

mw <- function(L, a, b) mw_matrix(L[a, , drop = FALSE], L[b, , drop = FALSE])

# log2 ~ группа + ковариаты (complete case), МНК для всех признаков одновременно
ols_group <- function(L, meta, case, ref, covariates) {
  idx <- c(case, ref)
  D <- data.frame(group = as.numeric(idx %in% case), row.names = idx)
  for (c in covariates) D[[c]] <- if (c == "sex") c(M = 1, F = 0)[meta[idx, "sex"]] else as.numeric(meta[idx, c])
  D <- D[stats::complete.cases(D), , drop = FALSE]
  used <- covariates[vapply(covariates, function(c) stats::sd(D[[c]]) > 0, TRUE)]
  Xd <- cbind(1, as.matrix(D[, c("group", used), drop = FALSE]))
  Y <- L[rownames(D), , drop = FALSE]
  XtX <- MASS_ginv(crossprod(Xd))
  beta <- XtX %*% crossprod(Xd, Y)
  res <- Y - Xd %*% beta
  dof <- nrow(Xd) - ncol(Xd)
  se <- sqrt(colSums(res^2) / dof * XtX[2, 2])
  p <- 2 * stats::pt(abs(beta[2, ] / se), dof, lower.tail = FALSE)
  list(beta = beta[2, ], p = p, n_case = sum(D$group), n_ref = sum(1 - D$group), used = used)
}

MASS_ginv <- function(A) {          # псевдообратная матрица (как numpy.linalg.pinv)
  s <- svd(A); tol <- max(dim(A)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# Оптимальное 1:1 назначение (венгерский алгоритм, clue::solve_LSAP) по |Δ| ≤ caliper внутри пола
optimal_pairs <- function(meta, case, ref, value, caliper, within_sex = TRUE) {
  strata <- if (within_sex) c("M", "F") else NA
  out <- list()
  for (sx in strata) {
    a <- case[!is.na(value[case]) & (is.na(sx) | meta[case, "sex"] == sx)]
    b <- ref[!is.na(value[ref]) & (is.na(sx) | meta[ref, "sex"] == sx)]
    if (!length(a) || !length(b)) next
    D <- abs(outer(value[a], value[b], "-"))
    cost <- ifelse(D <= caliper, D, 1e6)
    tr <- nrow(cost) > ncol(cost)
    asg <- as.integer(clue::solve_LSAP(if (tr) t(cost) else cost))
    ij <- if (tr) cbind(asg, seq_along(asg)) else cbind(seq_along(asg), asg)
    ok <- cost[ij] < 1e6
    if (any(ok)) out[[length(out) + 1]] <- data.frame(
      case = a[ij[ok, 1]], ref = b[ij[ok, 2]], sex = ifelse(is.na(sx), "", sx),
      value_case = value[a[ij[ok, 1]]], value_ref = value[b[ij[ok, 2]]], diff = D[ij[ok, , drop = FALSE]],
      stringsAsFactors = FALSE)
  }
  if (!length(out)) return(data.frame(case = character(0), ref = character(0), sex = character(0),
                                      value_case = numeric(0), value_ref = numeric(0), diff = numeric(0)))
  do.call(rbind, out)
}

age_match <- function(meta, case, ref, cfg)
  optimal_pairs(meta, case, ref, stats::setNames(meta$age, meta$sample), cfg$design$match_caliper_years,
                cfg$design$match_within_sex)

contam_match <- function(meta, case, ref, cfg, caliper_sd = 0.25) {
  z <- stats::setNames(meta$contam_index, meta$sample)
  optimal_pairs(meta, case, ref, z, caliper_sd * stats::sd(z[c(case, ref)], na.rm = TRUE), cfg$design$match_within_sex)
}

run_comparison <- function(L, meta, cm, comp_cat, cfg) {
  s <- comparison_sets(meta, cm); case <- s$case; ref <- s$ref
  out <- data.frame(row.names = colnames(L))
  info <- list(comparison = cm$name, n_case = length(case), n_ref = length(ref))
  r <- mw(L, case, ref)
  out$log2FC <- r$fc; out$AUC <- r$auc; out$p_MW <- r$p; out$q_MW <- bh(r$p)
  covs <- intersect(cfg$design$covariates, c(names(meta)))
  m <- ols_group(L, meta, case, ref, covs)
  out$beta_model <- m$beta; out$p_model <- m$p; out$q_model <- bh(m$p)
  info$model_covariates <- paste(m$used, collapse = "+"); info$model_n_case <- m$n_case; info$model_n_ref <- m$n_ref
  pairs <- age_match(meta, case, ref, cfg)
  r <- mw(L, pairs$case, pairs$ref)
  out$log2FC_matched <- r$fc; out$AUC_matched <- r$auc; out$p_matched <- r$p; out$q_matched <- bh(r$p)
  info$matched_pairs <- nrow(pairs)
  if (cm$ref == "Control" && !is.null(meta$contam_level)) {
    lo <- ref[meta[ref, "contam_level"] == "low"]; hi <- ref[meta[ref, "contam_level"] == "high"]
    r <- mw(L, case, lo); out$log2FC_lowindex <- r$fc; out$p_lowindex <- r$p
    r <- mw(L, hi, lo); out$log2FC_ctrl_hi_vs_lo <- r$fc; out$p_ctrl_hi_vs_lo <- r$p
    info$n_ref_lowindex <- length(lo); info$n_ref_highindex <- length(hi)
  }
  out$compound_category <- comp_cat[rownames(out)]
  out$finding <- finding_rule(out, cfg)
  info$n_q_model <- sum(out$q_model < cfg$stats$fdr, na.rm = TRUE)
  info$n_findings <- sum(out$finding)
  list(table = out, info = as.data.frame(info, stringsAsFactors = FALSE), pairs = cbind(comparison = cm$name, pairs))
}

# «Находка»: q(модель) < FDR и p(пары по возрасту) < confirm_p в том же направлении; для сравнений
# с контролем ещё p(против контролей с низким индексом) < confirm_p в том же направлении и
# отсутствие различия между контролями с высоким и низким индексом. Экспозиции находками не бывают.
finding_rule <- function(r, cfg) {
  s <- cfg$stats; sg <- sign(r$beta_model)
  ok <- r$q_model < s$fdr & r$p_matched < s$confirm_p & sign(r$log2FC_matched) == sg
  if (!is.null(r$p_lowindex))
    ok <- ok & r$p_lowindex < s$confirm_p & sign(r$log2FC_lowindex) == sg & r$p_ctrl_hi_vs_lo >= s$ctrl_index_p
  ok <- ok & !(r$compound_category %in% EXPOSURE)
  ok[is.na(ok)] <- FALSE
  ok
}

demographics <- function(meta, groups) {
  s <- meta[meta$study, ]
  tab <- do.call(rbind, lapply(groups, function(g) {
    d <- s[s$group == g, ]
    data.frame(group = g, n = nrow(d), male = sum(d$sex == "M"), female = sum(d$sex == "F"),
               age_n = sum(!is.na(d$age)), age_mean = mean(d$age, na.rm = TRUE), age_sd = stats::sd(d$age, na.rm = TRUE),
               age_median = stats::median(d$age, na.rm = TRUE), age_min = suppressWarnings(min(d$age, na.rm = TRUE)),
               age_max = suppressWarnings(max(d$age, na.rm = TRUE)))
  }))
  tests <- list()
  ok_age <- groups[vapply(groups, function(g) sum(!is.na(s$age[s$group == g])) > 1, TRUE)]
  if (length(ok_age) > 1)
    tests[[1]] <- data.frame(test = "возраст, Краскел–Уоллис", p = stats::kruskal.test(age ~ group, data = s[s$group %in% ok_age, ])$p.value)
  for (pr in utils::combn(groups, 2, simplify = FALSE)) {
    x <- s$age[s$group == pr[1]]; y <- s$age[s$group == pr[2]]
    if (sum(!is.na(x)) > 1 && sum(!is.na(y)) > 1)
      tests[[length(tests) + 1]] <- data.frame(test = sprintf("возраст %s vs %s, Манн–Уитни", pr[1], pr[2]),
        p = stats::wilcox.test(x, y, exact = FALSE, correct = TRUE)$p.value)
    d <- s[s$group %in% pr & s$sex %in% c("M", "F"), ]
    ct <- table(factor(d$group, levels = pr), factor(d$sex, levels = c("F", "M")))
    tests[[length(tests) + 1]] <- data.frame(test = sprintf("пол %s vs %s, точный тест Фишера", pr[1], pr[2]),
                                             p = stats::fisher.test(ct)$p.value)
  }
  list(table = tab, tests = do.call(rbind, tests))
}

# Модули 2–4: контроль качества, фильтр устойчивых признаков, подготовка матрицы образцов.

rsd_cols <- function(M) apply(M, 2, function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))
half_min_log2 <- function(X) { X <- apply(X, 2, function(v) { v[is.na(v)] <- min(v, na.rm = TRUE) / 2; v }); log2(X) }

pca_plot <- function(L, groups, title, ellipses = TRUE) {
  L <- L[, apply(L, 2, stats::sd) > 0, drop = FALSE]
  p <- stats::prcomp(L, center = TRUE, scale. = TRUE)
  v <- 100 * p$sdev^2 / sum(p$sdev^2)
  d <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], group = groups)
  g <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_colour_manual(values = group_colors(d$group)) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", v[1]), y = sprintf("PC2 (%.1f%%)", v[2]), colour = NULL, title = title) + theme_um()
  if (ellipses) g <- g + ggplot2::stat_ellipse(data = d[d$group != "QC" & ave(seq_along(d$group), d$group, FUN = length) > 3, ],
                                               level = 0.95, linewidth = 0.3)
  g
}

# ------------------------------------------------------------------ модуль 2
module2_qc <- function(s) {
  out <- module_dir(s, 2); d <- load_previous(s, 2)
  X <- d$X; meta <- d$meta[rownames(X), ]
  smp <- meta$kind == "sample"; qc <- meta$kind == "QC"
  summ <- data.frame(parameter = c("Инъекций всего", "из них QC", "Признаков", "Пропусков, все инъекции, %",
                                   "Пропусков в образцах, %", "Пропусков в QC, %", "Медианный RSD в QC, %"),
                     value = c(nrow(X), sum(qc), ncol(X), round(100 * mean(is.na(X)), 1),
                               round(100 * mean(is.na(X[smp, ])), 1), if (any(qc)) round(100 * mean(is.na(X[qc, ])), 1) else NA,
                               if (sum(qc) >= 2) round(100 * stats::median(rsd_cols(X[qc, , drop = FALSE]), na.rm = TRUE), 1) else NA))
  by_group <- stats::aggregate(list(missing_pct = rowMeans(is.na(X)) * 100), list(group = meta$group), mean)
  write_csv(summ, file.path(out, "qc_summary.csv")); write_csv(by_group, file.path(out, "missing_by_group.csv"))
  for (i in seq_len(nrow(summ))) log_msg("  ", summ$parameter[i], ": ", summ$value[i])
  ok <- colSums(!is.na(X)) >= 2
  save_fig(s, pca_plot(half_min_log2(X[, ok]), meta$group, "PCA: образцы и QC (log2, автошкалирование)"), out, "pca_qc", 12, 11)
  tic <- data.frame(order = meta$injection_order, tic = rowSums(X, na.rm = TRUE), group = meta$group)
  g <- ggplot2::ggplot(tic, ggplot2::aes(order, tic, colour = group, shape = group == "QC")) + ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_y_log10() + ggplot2::scale_colour_manual(values = group_colors(tic$group)) + ggplot2::guides(shape = "none") +
    ggplot2::labs(x = "Порядок ввода", y = "Суммарная интенсивность (TIC)", colour = NULL) + theme_um()
  save_fig(s, g, out, "tic_by_order", 14, 9)
  if (sum(qc) >= 2) {
    r <- data.frame(rsd = 100 * rsd_cols(X[qc, , drop = FALSE])); r <- r[is.finite(r$rsd), , drop = FALSE]
    g <- ggplot2::ggplot(r, ggplot2::aes(rsd)) + ggplot2::geom_histogram(bins = 50, fill = "grey60", colour = "white") +
      ggplot2::geom_vline(xintercept = 100 * s$filter$max_qc_rsd, linetype = 2) +
      ggplot2::labs(x = "RSD в QC, %", y = "Число признаков") + theme_um()
    save_fig(s, g, out, "rsd_qc", 10, 8)
  }
  save_module(s, 2, c(d, list(qc_summary = summ)))
}

# ------------------------------------------------------------------ модуль 3
# Два профиля фильтра. Стандартный (по умолчанию для анализа, модули 4–8): присутствие ≥ 50 % хотя бы
# в одной группе, обнаружение в QC ≥ 80 %, бланк, RSD в QC ≤ 30 %. Строгий (для базы признаков, модуль 9):
# обнаружение во всех образцах, в QC ≥ 80 %, бланк, RSD в QC ≤ 15 % — сигналы, устойчивые в любой когорте
# и пригодные для накопления; признаки, характерные только для одной группы, в него не входят.
module3_filter <- function(s) {
  out <- module_dir(s, 3); d <- load_module(s, 2); f <- s$filter
  X <- d$X; meta <- d$meta[rownames(X), ]
  qc <- meta$kind == "QC"; bl <- meta$kind == "blank"; smp <- meta$kind == "sample"
  groups <- sort(unique(meta$group[smp]))
  tab <- data.frame(id = colnames(X), mz = d$info$mz, rt_min = d$info$rt_min, stringsAsFactors = FALSE)
  for (g in groups) tab[[paste0("detect_", g)]] <- colMeans(!is.na(X[meta$group == g, , drop = FALSE]))
  tab$detect_all <- colMeans(!is.na(X[smp, , drop = FALSE]))
  if (sum(qc) >= 3) { tab$detect_QC <- colMeans(!is.na(X[qc, , drop = FALSE])); tab$rsd_qc <- rsd_cols(X[qc, , drop = FALSE]) }
  else log_msg("  QC меньше трёх — фильтры по QC не применяются")
  if (any(bl)) {
    B <- X[bl, , drop = FALSE]; B[is.na(B)] <- 0
    ref <- if (any(qc)) X[qc, , drop = FALSE] else X[smp, , drop = FALSE]
    mb <- apply(B, 2, stats::median); mb[mb == 0] <- NA
    tab$blank_ratio <- apply(ref, 2, stats::median, na.rm = TRUE) / mb
    tab$blank_ratio[is.na(tab$blank_ratio)] <- Inf
  }
  run_profile <- function(name, presence, presence_label, min_qc, max_rsd) {
    st <- data.frame(profile = name, step = "Исходно", n = ncol(X)); keep <- presence
    add <- function(lbl) st[nrow(st) + 1, ] <<- list(name, lbl, sum(keep))
    add(presence_label)
    if (sum(qc) >= 3) { keep <- keep & tab$detect_QC >= min_qc; add(sprintf("Обнаружение в QC ≥ %.0f%%", 100 * min_qc)) }
    if (any(bl)) { keep <- keep & tab$blank_ratio >= f$blank_ratio; add(sprintf("Отношение к бланку ≥ %g", f$blank_ratio)) }
    if (sum(qc) >= 3) { keep <- keep & !is.na(tab$rsd_qc) & tab$rsd_qc <= max_rsd; add(sprintf("RSD в QC ≤ %.0f%%", 100 * max_rsd)) }
    list(keep = keep, steps = st)
  }
  ps <- run_profile("стандартный", apply(tab[paste0("detect_", groups)], 1, max) >= f$min_group_detect,
                    sprintf("Присутствие ≥ %.0f%% хотя бы в одной группе", 100 * f$min_group_detect), f$min_qc_detect, f$max_qc_rsd)
  sp <- f$strict
  pt <- run_profile("строгий", tab$detect_all >= sp$min_all_detect,
                    sprintf("Обнаружение в ≥ %.0f%% всех образцов", 100 * sp$min_all_detect), sp$min_qc_detect, sp$max_qc_rsd)
  tab$kept_standard <- ps$keep; tab$kept_strict <- pt$keep
  steps <- rbind(ps$steps, pt$steps)
  for (i in seq_len(nrow(steps))) log_msg(sprintf("  [%s] %s: %d", steps$profile[i], steps$step[i], steps$n[i]))
  inr <- function(n) n >= f$target_min && n <= f$target_max
  log_msg(sprintf("  Устойчивых признаков: стандартный профиль %d, строгий %d (целевой диапазон базы признаков %d–%d: %s)",
                  sum(ps$keep), sum(pt$keep), f$target_min, f$target_max, if (inr(sum(pt$keep))) "строгий профиль в пределах" else "строгий профиль вне"))
  keep <- if (f$analysis_profile == "strict") pt$keep else ps$keep
  log_msg(sprintf("  Для анализа (модули 4–8) используется %s профиль: %d признаков; для базы признаков (модуль 9) — %s",
                  if (f$analysis_profile == "strict") "строгий" else "стандартный", sum(keep),
                  if (s$report$database_profile == "strict") "строгий" else "стандартный"))
  if (sum(keep) < 2) stop("После фильтрации осталось меньше двух признаков — ослабьте пороги в настройках")
  write_csv(tab, file.path(out, "feature_filter.csv")); write_csv(steps, file.path(out, "filter_steps.csv"))
  Xf <- X[, keep, drop = FALSE]
  write_csv(cbind(file = rownames(Xf), as.data.frame(Xf)), file.path(out, "filtered_peak_table.csv"))
  Xs <- X[, pt$keep, drop = FALSE]
  write_csv(cbind(file = rownames(Xs), as.data.frame(Xs)), file.path(out, "strict_profile_peak_table.csv"))
  save_module(s, 3, list(X = Xf, info = d$info[keep, ], meta = meta, filter = tab, steps = steps,
                         strict_ids = tab$id[pt$keep], standard_ids = tab$id[ps$keep], in_target = inr(sum(pt$keep))))
}

# ------------------------------------------------------------------ модуль 4
module4_samples <- function(s) {
  out <- module_dir(s, 4); d <- load_module(s, 3)
  X <- d$X; meta <- d$meta[rownames(X), ]
  use <- meta$kind == "sample"
  if (length(s$samples$groups)) use <- use & meta$group %in% unlist(s$samples$groups)
  X <- X[use, , drop = FALSE]; meta <- meta[use, ]
  log_msg(sprintf("  Удалены QC и бланки: осталось %d инъекций", nrow(X)))
  if (isTRUE(s$samples$average_replicates) && any(duplicated(meta$sample))) {
    ng <- tapply(meta$group, meta$sample, function(g) length(unique(g)))
    if (any(ng > 1)) stop("У повторов одного образца разные группы: ", paste(names(ng)[ng > 1], collapse = ", "))
    sm <- unique(meta$sample)
    A <- t(vapply(sm, function(k) { v <- colMeans(X[meta$sample == k, , drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA; v }, numeric(ncol(X))))
    m2 <- meta[match(sm, meta$sample), ]
    m2$injection_order <- tapply(meta$injection_order, meta$sample, min)[sm]
    m2$n_injections <- as.integer(table(meta$sample)[sm])
    log_msg(sprintf("  Усреднены повторы: %d инъекций → %d образцов", nrow(X), length(sm)))
    X <- A; rownames(X) <- sm; meta <- m2; rownames(meta) <- sm
  } else { rownames(X) <- meta$sample; rownames(meta) <- meta$sample; meta$n_injections <- 1L }
  o <- order(meta$injection_order)
  X <- X[o, , drop = FALSE]; meta <- meta[o, ]
  tab <- data.frame(sample = rownames(X), group = meta$group, injection_order = meta$injection_order, as.data.frame(X), check.names = FALSE)
  write_csv(tab, file.path(out, "samples_matrix.csv"))
  log_msg(sprintf("  Матрица образцов: %d × %d; группы: %s", nrow(X), ncol(X), paste(sprintf("%s %d", names(table(meta$group)), table(meta$group)), collapse = ", ")))
  save_module(s, 4, list(X = X, info = d$info, meta = meta))
}

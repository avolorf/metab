# Модуль 8. Визуализация: PCA и иерархическая кластеризация (Манхэттен, Уорд) на трёх этапах обработки,
# тепловая карта, диаграммы размаха по метаболитам, PLS-DA с VIP. PNG/JPEG с разрешением из настроек.

feature_label <- function(info, ids) {
  i <- match(ids, info$id)
  ifelse(nzchar(info$name[i]) & info$name[i] != "Unknown", substr(info$name[i], 1, 28),
         ifelse(is.na(info$mz[i]), ids, sprintf("m/z %.4f @ %.2f min", info$mz[i], info$rt_min[i])))
}

draw_dendrogram <- function(L, groups, title, s) {
  hc <- stats::hclust(stats::dist(L, method = s$visual$hca_distance), method = s$visual$hca_method)
  cols <- group_colors(groups); gi <- stats::setNames(as.character(groups), rownames(L))
  dd <- stats::dendrapply(stats::as.dendrogram(hc), function(n) {
    if (stats::is.leaf(n)) attr(n, "nodePar") <- list(lab.col = cols[gi[attr(n, "label")]], pch = NA, lab.cex = 0.35)
    n })
  graphics::par(mar = c(3, 3, 2, 0.5), cex.main = 0.9)
  plot(dd, main = title, ylab = "")
  k <- length(unique(groups))
  stats::rect.hclust(hc, k = k, border = "grey40")
  graphics::legend("topright", legend = names(cols), text.col = cols, bty = "n", cex = 0.7)
}

vip_pls <- function(m) {
  W <- m$loading.weights; Tt <- m$scores; Q <- m$Yloadings
  SS <- colSums(Tt^2) * colSums(matrix(Q^2, ncol = ncol(W)))
  Wn <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  drop(sqrt(nrow(W) * (Wn^2 %*% SS) / sum(SS)))
}

module8_visual <- function(s) {
  out <- module_dir(s, 8); d <- load_module(s, 7); v <- s$visual
  meta <- d$meta; g <- meta$group; info <- d$info; sel <- d$sel
  Ls <- d$L_sel %||% d$L                                  # данные, на которых выполнен отбор
  sets <- list("A — до коррекции" = d$L_raw, "B — после коррекции" = d$L,
               "C — после отбора (ОСФ+СПФ+РОП)" = Ls[, sel$f_rfe, drop = FALSE])
  pcs <- lapply(names(sets), function(n) if (ncol(sets[[n]]) >= 2) pca_plot(sets[[n]], g, n) else
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "менее двух признаков") + ggplot2::theme_void() + ggplot2::labs(title = n))
  save_fig(s, patchwork::wrap_plots(pcs, nrow = 1), out, "pca_stages", 26, 9)
  save_base_fig(s, function() { graphics::par(mfrow = c(1, 3))
    for (n in names(sets)) if (ncol(sets[[n]]) >= 1) draw_dendrogram(sets[[n]], g, n, s) }, out, "hca_stages", 30, 11)
  hf <- if (length(sel$f_spf) >= 2) sel$f_spf else utils::head(sel$uvf$id[order(sel$uvf$p)], v$heatmap_top)
  hf <- utils::head(hf, v$heatmap_top)
  M <- t(scale(Ls[, hf, drop = FALSE])); rownames(M) <- make.unique(feature_label(info, hf)); M[M > 3] <- 3; M[M < -3] <- -3
  ann <- data.frame(group = g, row.names = rownames(Ls))
  save_base_fig(s, function() pheatmap::pheatmap(M, annotation_col = ann, show_colnames = FALSE, fontsize_row = 5,
    clustering_distance_cols = v$hca_distance, clustering_method = v$hca_method, annotation_colors = list(group = group_colors(g)),
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), main = "Heatmap: z-score, Ward"), out, "heatmap", 18, 16)
  bf <- unique(c(sel$f_rfe, utils::head(sel$uvf$id[order(sel$uvf$p)], v$boxplot_top)))[seq_len(min(v$boxplot_top, ncol(Ls)))]
  bd <- do.call(rbind, lapply(bf, function(f) data.frame(name = feature_label(info, f), value = Ls[, f], group = g)))
  bd$name <- factor(bd$name, levels = unique(bd$name))
  bp <- ggplot2::ggplot(bd, ggplot2::aes(group, value, fill = group)) + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 0.5) + ggplot2::facet_wrap(~name, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = group_colors(g), guide = "none") + ggplot2::labs(x = NULL, y = "log2 интенсивность") + theme_um() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6))
  save_fig(s, bp, out, "boxplots_metabolites", 17, 4.5 * ceiling(length(bf) / 4) + 1)
  Y <- stats::model.matrix(~ factor(g) - 1); Xs <- scale(d$L)
  m <- pls::plsr(Y ~ Xs, ncomp = 2, scale = FALSE, method = if (ncol(Xs) > nrow(Xs)) "widekernelpls" else "kernelpls")
  vip <- vip_pls(m); vt <- utils::head(order(-vip), v$vip_top)
  vd <- data.frame(name = factor(feature_label(info, colnames(d$L)[vt]), levels = rev(make.unique(feature_label(info, colnames(d$L)[vt])))), vip = vip[vt])
  vd$name <- factor(make.unique(as.character(vd$name)), levels = rev(make.unique(as.character(vd$name))))
  pv <- ggplot2::ggplot(vd, ggplot2::aes(vip, name)) + ggplot2::geom_col(fill = "#0072B2") + ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::labs(x = "VIP", y = NULL, title = "PLS-DA: важность переменных (VIP)") + theme_um()
  sc <- data.frame(t1 = m$scores[, 1], t2 = m$scores[, 2], group = g)
  ps <- ggplot2::ggplot(sc, ggplot2::aes(t1, t2, colour = group)) + ggplot2::geom_point(size = 1.6) + ggplot2::stat_ellipse(level = 0.95) +
    ggplot2::scale_colour_manual(values = group_colors(g)) + ggplot2::labs(x = "Компонента 1", y = "Компонента 2", colour = NULL, title = "PLS-DA: счета (описательно)") + theme_um()
  save_fig(s, patchwork::wrap_plots(ps, pv, nrow = 1, widths = c(1, 1.2)), out, "plsda_vip", 24, 10)
  write_csv(merge(info[, c("id", "mz", "rt_min", "name")], data.frame(id = colnames(d$L), vip = vip), by = "id")[order(-vip), ], file.path(out, "plsda_vip.csv"))
  log_msg(sprintf("  Рисунки: PCA и ИКА по трём этапам, тепловая карта (%d признаков), диаграммы размаха (%d), PLS-DA/VIP; форматы: %s, %d dpi",
                  length(hf), length(bf), paste(v$formats, collapse = ", "), v$dpi))
  save_module(s, 8, list(figures = list.files(out, pattern = "\\.(png|jpg)$")))
}

# Визуализация: контроль качества, распределения интенсивностей, PCA, индекс загрязнения,
# вулканы, тепловые карты, диаграммы размаха, ROC.

GROUP_COLORS <- c(Control = "#0072B2", BC = "#D55E00", PC = "#009E73", QC = "#7F7F7F", Blank = "#CC79A7")
group_cols <- function(g) {
  u <- unique(as.character(g)); col <- GROUP_COLORS[u]
  extra <- is.na(col)
  if (any(extra)) col[extra] <- c("#E69F00", "#56B4E9", "#F0E442", "#000000")[seq_len(sum(extra))]
  stats::setNames(col, u)
}

save_plot <- function(p, path, w = 17, h = 10, dpi = 300) {
  ggplot2::ggsave(paste0(path, ".png"), p, width = w, height = h, units = "cm", dpi = dpi, bg = "white")
  invisible(paste0(path, ".png"))
}

theme_um <- function() ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top")

plot_qc_overview <- function(X, meta, rsd_qc, factors, path, dpi) {
  m <- meta[!is.na(meta$inj_order) & (meta$study | meta$kind == "QC"), ]
  d <- data.frame(order = m$inj_order, tic = rowSums(X[m$sample, ], na.rm = TRUE), group = m$group)
  a <- ggplot2::ggplot(d, ggplot2::aes(order, tic, colour = group, shape = group == "QC")) +
    ggplot2::geom_point(size = 1.6) + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = group_cols(d$group)) + ggplot2::guides(shape = "none") +
    ggplot2::labs(x = "Порядок инъекции", y = "Суммарная интенсивность", colour = NULL, title = "(a) TIC по порядку ввода") + theme_um()
  b <- ggplot2::ggplot(data.frame(rsd = 100 * rsd_qc), ggplot2::aes(rsd)) +
    ggplot2::geom_histogram(bins = 50, fill = "grey60", colour = "white") +
    ggplot2::geom_vline(xintercept = 30, linetype = 2) +
    ggplot2::labs(x = "RSD в QC, %", y = "Число признаков", title = "(b) RSD признаков в QC") + theme_um()
  f <- data.frame(PQN = factors$PQN, MSTUS = factors$MSTUS, group = meta[rownames(factors), "group"])
  c <- ggplot2::ggplot(f, ggplot2::aes(PQN, MSTUS, colour = group)) + ggplot2::geom_point(size = 1.4) +
    ggplot2::geom_abline(linetype = 3) + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = group_cols(f$group)) +
    ggplot2::labs(title = "(c) Факторы разведения", colour = NULL) + theme_um()
  save_plot(patchwork::wrap_plots(a, b, c, nrow = 1, widths = c(1.4, 1, 1)), path, 24, 8, dpi)
}

# Диаграммы размаха log2-интенсивностей по пробам (в порядке ввода) до и после нормализации
plot_sample_distributions <- function(Lraw, Lnorm, meta, path, dpi) {
  one <- function(L, title) {
    ids <- rownames(L)[order(meta[rownames(L), "inj_order"])]
    d <- data.frame(sample = factor(rep(ids, ncol(L)), levels = ids), value = as.vector(L[ids, ]),
                    group = rep(meta[ids, "group"], ncol(L)))
    d <- d[is.finite(d$value), ]
    ggplot2::ggplot(d, ggplot2::aes(sample, value, fill = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = group_cols(d$group)) +
      ggplot2::labs(x = "Пробы (в порядке ввода)", y = "log2 интенсивность", fill = NULL, title = title) +
      theme_um() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  }
  save_plot(patchwork::wrap_plots(one(Lraw, "(a) до нормализации"), one(Lnorm, "(b) после нормализации и импутации"), ncol = 1),
            path, 24, 14, dpi)
}

plot_pca <- function(pc, meta, path, dpi) {
  d <- data.frame(pc$scores[, 1:2], group = meta[rownames(pc$scores), "group"],
                  index = if (!is.null(meta$contam_index)) meta[rownames(pc$scores), "contam_index"] else NA)
  lab <- sprintf("PC%d (%.1f%%)", 1:2, 100 * pc$var[1:2])
  a <- ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, colour = group)) + ggplot2::geom_point(size = 1.6) +
    ggplot2::stat_ellipse(data = d[d$group != "QC", ], level = 0.95, linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = group_cols(d$group)) +
    ggplot2::labs(x = lab[1], y = lab[2], colour = NULL, title = "(a) PCA, Pareto, группы и QC") + theme_um()
  ps <- list(a)
  if (!all(is.na(d$index))) ps[[2]] <- ggplot2::ggplot(d[d$group != "QC", ], ggplot2::aes(PC1, PC2, colour = index)) +
    ggplot2::geom_point(size = 1.6) + ggplot2::scale_colour_viridis_c() +
    ggplot2::labs(x = lab[1], y = lab[2], colour = "индекс", title = "(b) индекс загрязнения") + theme_um()
  save_plot(patchwork::wrap_plots(ps, nrow = 1), path, 8.5 * length(ps), 8.5, dpi)
}

plot_contamination <- function(meta, groups, path, dpi) {
  d <- meta[meta$study, ]; d$group <- factor(d$group, levels = groups)
  p <- ggplot2::ggplot(d, ggplot2::aes(group, contam_index, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) + ggplot2::geom_jitter(width = 0.15, size = 1) +
    ggplot2::scale_fill_manual(values = group_cols(d$group), guide = "none") +
    ggplot2::labs(x = NULL, y = "Индекс загрязнения (средний z)") + theme_um()
  save_plot(p, path, 8.5, 8, dpi)
}

plot_volcano <- function(r, title, names, fdr, path, dpi) {
  d <- data.frame(fc = r$beta_model, p = -log10(r$p_model), q = r$q_model, find = r$finding,
                  expo = r$compound_category %in% EXPOSURE, name = names[rownames(r)])
  d$cls <- ifelse(d$expo, "экспозиция", ifelse(d$find, "находка", ifelse(d$q < fdr, "q < FDR", "н/з")))
  thr <- if (any(d$q < fdr, na.rm = TRUE)) min(d$p[d$q < fdr], na.rm = TRUE) else NA
  lab <- d[d$find, ]; lab <- utils::head(lab[order(-lab$p), ], 10)
  p <- ggplot2::ggplot(d, ggplot2::aes(fc, p, colour = cls, shape = cls)) + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = c("н/з" = "grey75", "q < FDR" = "#56B4E9", "находка" = "#D55E00", "экспозиция" = "black")) +
    ggplot2::scale_shape_manual(values = c("н/з" = 16, "q < FDR" = 16, "находка" = 16, "экспозиция" = 4)) +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = name), size = 2.2, colour = "black", max.overlaps = 20) +
    ggplot2::labs(x = "log2FC (модель с ковариатами)", y = "-log10 p", colour = NULL, shape = NULL, title = title) + theme_um()
  if (!is.na(thr)) p <- p + ggplot2::geom_hline(yintercept = thr, linetype = 2)
  save_plot(p, path, 12, 10, dpi)
}

plot_heatmap <- function(L, meta, r, case, ref, names, top, path) {
  s <- c(meta$sample[meta$study & meta$group == case], meta$sample[meta$study & meta$group == ref])
  ok <- !(r$compound_category %in% EXPOSURE)
  f <- rownames(r)[ok][order(r$p_model[ok])][seq_len(min(top, sum(ok)))]
  if (length(f) < 2) return(invisible(NULL))
  M <- t(scale(L[s, f, drop = FALSE])); rownames(M) <- names[f]
  M[M > 3] <- 3; M[M < -3] <- -3
  ann <- data.frame(group = meta[s, "group"], row.names = s)
  grDevices::png(paste0(path, ".png"), width = 18, height = 16, units = "cm", res = 300,
                 type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
  pheatmap::pheatmap(M, annotation_col = ann, show_colnames = FALSE, fontsize_row = 5,
                     clustering_distance_rows = "correlation", clustering_method = "ward.D2",
                     color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
                     annotation_colors = list(group = group_cols(ann$group)),
                     main = sprintf("top-%d compounds by p (z-score, clipped at 3)", length(f)))  # pheatmap: ASCII-заголовок
  grDevices::dev.off()
  invisible(paste0(path, ".png"))
}

plot_boxplots <- function(L, meta, feats, names, groups, title, path, dpi) {
  if (!length(feats)) return(invisible(NULL))
  s <- meta$sample[meta$study]
  d <- do.call(rbind, lapply(feats, function(f) data.frame(name = unname(names[f]), value = unname(L[s, f]), group = meta[s, "group"])))
  d$group <- factor(d$group, levels = groups); d$name <- factor(d$name, levels = unique(names[feats]))
  p <- ggplot2::ggplot(d, ggplot2::aes(group, value, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) + ggplot2::geom_jitter(width = 0.15, size = 0.5) +
    ggplot2::facet_wrap(~name, scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = group_cols(d$group), guide = "none") +
    ggplot2::labs(x = NULL, y = "log2 интенсивность (норм.)", title = title) + theme_um() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 6))
  save_plot(p, path, 17, 4.5 * ceiling(length(feats) / 4) + 1.5, dpi)
}

plot_roc <- function(curves, title, path, dpi) {
  if (!length(curves)) return(invisible(NULL))
  d <- do.call(rbind, lapply(curves, function(c) data.frame(fpr = FPR_GRID, tpr = c$tpr,
                                                          model = sprintf("%s, AUC %.2f", c$label, c$auc))))
  p <- ggplot2::ggplot(d, ggplot2::aes(fpr, tpr, colour = model)) + ggplot2::geom_line() +
    ggplot2::geom_abline(linetype = 3) + ggplot2::coord_equal() +
    ggplot2::labs(x = "1 − специфичность", y = "Чувствительность", colour = NULL, title = title) +
    theme_um() + ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 6))
  save_plot(p, path, 20, 9, dpi)
}

plot_norm_concordance <- function(a, b, title, path, dpi) {
  d <- data.frame(x = a$beta_model, y = b$beta_model)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point(size = 0.5, alpha = 0.4) + ggplot2::geom_abline(linetype = 3) +
    ggplot2::labs(x = "log2FC, PQN", y = "log2FC, MSTUS",
                  title = sprintf("%s: ρ = %.2f", title, stats::cor(d$x, d$y, method = "spearman", use = "complete.obs"))) + theme_um()
  save_plot(p, path, 8.5, 8.5, dpi)
}

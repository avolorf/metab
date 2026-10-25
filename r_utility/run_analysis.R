#!/usr/bin/env Rscript
# UroMetabR — полуавтоматическая обработка данных ненаправленного ВЭЖХ-МС профилирования мочи.
#
#   Rscript run_analysis.R config/config.yaml           # полный прогон
#   Rscript run_analysis.R config/config.yaml --quick   # быстрая проверка (3 повтора CV, 20 перестановок)
#
# Из RStudio: открыть файл, задать путь к config.yaml в строке cfg_path ниже и нажать Source.

main <- function(cfg_path = NULL, quick = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  if (is.null(cfg_path) && length(args)) cfg_path <- args[!startsWith(args, "--")][1]
  quick <- quick || "--quick" %in% args
  here <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))),
                   error = function(e) getwd())
  if (is.na(here) || !dir.exists(file.path(here, "R"))) here <- getwd()
  if (!isTRUE(l10n_info()$`UTF-8`))                     # модули и подписи в UTF-8
    for (loc in c("C.UTF-8", "en_US.UTF-8", "ru_RU.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  options(urometabr.dir = file.path(here, "R"))
  for (f in sort(list.files(file.path(here, "R"), pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
  check_packages()
  cfg <- read_config(cfg_path, quick)
  run_pipeline(cfg)
}

check_packages <- function() {
  need <- c("data.table", "yaml", "ggplot2", "patchwork", "ggrepel", "pheatmap", "openxlsx", "pls", "clue",
            "jsonlite", "digest")
  miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss)) stop("Не установлены пакеты: ", paste(miss, collapse = ", "),
                         "\nУстановка: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")
}

run_pipeline <- function(cfg) {
  t0 <- Sys.time(); log_reset()
  set.seed(cfg$seed)
  out <- cfg$output
  if (dir.exists(out)) unlink(out, recursive = TRUE)
  dir.create(file.path(out, "figures"), recursive = TRUE)
  fig <- function(n) file.path(out, "figures", n)
  dpi <- cfg$report$fig_dpi; groups <- cfg$design$groups; fdr <- cfg$stats$fdr
  sheets <- list(); html <- list()

  # ---------------------------------------------------------------- 1. данные
  log_msg("UroMetabR ", UROMETABR_VERSION, "; R ", R.version$major, ".", R.version$minor)
  ft <- read_feature_table(cfg)
  meta <- load_metadata(cfg$input$metadata, ft$kind, cfg)
  meta$inj_order <- injection_order(meta, ft$inj_order, cfg$design$qc_every)
  inj <- stats::setNames(meta$inj_order, meta$sample)
  study <- meta$sample[meta$study]; qc <- meta$sample[meta$kind == "QC"]
  if (length(qc) < 3) stop("Нужно не менее 3 инъекций QC (найдено ", length(qc), ")")
  log_msg(sprintf("Таблица: %d инъекций (%d QC, %d бланков), %d признаков; дубликатов строк удалено: %d",
                  nrow(ft$X), length(qc), sum(meta$kind == "blank"), ncol(ft$X), ft$n_dup))
  log_msg("Включено проб: ", paste(sprintf("%s %d", groups, vapply(groups, function(g) sum(meta$group[meta$study] == g), 0L)), collapse = ", "))
  log_msg(sprintf("Доля пропусков в исходной таблице (исследуемые пробы): %.1f%%", 100 * mean(is.na(ft$X[study, ]))))
  dm <- demographics(meta, groups)
  sheets$demographics <- dm$table; sheets$demographics_tests <- dm$tests

  # ---------------------------------------------------------------- 2. предобработка
  use <- meta$sample[meta$study | meta$kind %in% c("QC", "blank")]
  flt <- filter_features(ft$X[use, , drop = FALSE], meta[use, ], cfg)
  X1 <- flt$X; feat <- flt$table
  loo <- if (length(qc) >= 5) drift_loo(X1, inj, qc) else c(none = NA, lowess = NA)
  drift <- cfg$preprocess$drift
  if (drift == "auto") drift <- if (!is.na(loo["lowess"]) && loo["lowess"] < 0.8 * loo["none"]) "lowess" else "none"
  log_msg(sprintf("Дрейф (leave-one-QC-out, медиана |Δlog2|): без коррекции %.3f, LOWESS %.3f → коррекция: %s",
                  loo["none"], loo["lowess"], drift))
  if (drift == "lowess") X1 <- drift_correct(X1, inj, qc)
  ana <- meta$sample[meta$study | meta$kind == "QC"]
  norm <- list(); facs <- data.frame(row.names = ana)
  for (m in unique(c(cfg$preprocess$normalization, cfg$preprocess$compare_normalizations))) {
    n <- normalize(X1[ana, , drop = FALSE], meta, m)
    norm[[m]] <- impute_log(n$X); facs[[m]] <- n$f
    log_msg(sprintf("Нормализация %s: факторы проб %.2f–%.2f", n$info, min(n$f[study]), max(n$f[study])))
  }
  L <- norm[[cfg$preprocess$normalization]]
  meta$dilution <- NA; meta[ana, "dilution"] <- facs[[cfg$preprocess$normalization]]
  log_msg(sprintf("Основная нормализация: %s; импутировано %.1f%% значений (½ минимума), log2",
                  cfg$preprocess$normalization, 100 * mean(is.na(X1[ana, ]))))
  low <- study[meta[study, "dilution"] < cfg$preprocess$low_dilution_factor]
  if (length(low)) log_msg("Пробы с фактором разведения < ", cfg$preprocess$low_dilution_factor, " (проверить хроматограммы): ",
                           paste(low, collapse = ", "))
  for (g in groups) feat[[paste0("imputed_", g)]] <- 1 - feat[[paste0("det_", g)]]
  pa_rng <- cfg$preprocess$presence_absence
  feat$presence_absence <- presence_absence(feat, groups, pa_rng[1], pa_rng[2])

  # ---------------------------------------------------------------- 3. соединения, загрязнение, категории
  info <- ft$info[colnames(L), ]
  grp <- group_features(L[study, ], info, col_medians(X1[study, , drop = FALSE]), cfg$grouping$rt_tol, cfg$grouping$min_corr)
  reps <- grp$feature[grp$is_representative]
  log_msg(sprintf("Группировка ионов (ΔRT ≤ %g мин, ρ ≥ %g): %d признаков → %d соединений; многоионных %d",
                  cfg$grouping$rt_tol, cfg$grouping$min_corr, ncol(L), length(reps), sum(table(grp$compound) > 1)))
  series <- match_series(info, cfg$input$contaminant_series, cfg$grouping$ppm)
  ci <- if (!is.null(series)) contamination_index(L, series, study) else NULL
  if (!is.null(ci)) {
    meta$contam_index <- NA; meta[names(ci$index), "contam_index"] <- ci$index
    meta$contam_level <- NA; meta[names(ci$level), "contam_level"] <- ci$level
    log_msg(sprintf("Индекс загрязнения: %d членов ряда; средние: %s", length(ci$features),
                    paste(sprintf("%s %+.2f", groups, vapply(groups, function(g) mean(meta$contam_index[meta$study & meta$group == g]), 0)), collapse = ", ")))
  } else {
    log_msg("Ряд экзогенных соединений не задан или не найден: индекс загрязнения не используется")
    cfg$design$covariates <- setdiff(cfg$design$covariates, "contam_index")
  }
  cats <- categorize(info, grp, series, exposure_category(info, cfg$input$exposure_keywords))
  names_ <- display_name(info)
  comp_cat <- stats::setNames(cats[reps, "compound_category"], reps)
  feat <- cbind(feat[colnames(L), ], info[colnames(L), -1], grp[colnames(L), c("compound", "n_ions", "is_representative", "ion_relations")],
                cats[colnames(L), c("category", "compound_category")])
  log_msg("Категории соединений: ", paste(sprintf("%s %d", names(table(comp_cat)), as.integer(table(comp_cat))), collapse = ", "))

  # ---------------------------------------------------------------- 4. одномерная статистика
  Lc <- L[study, reps, drop = FALSE]
  res <- list(); summ <- list(); pairs_all <- list()
  for (cm in cfg$design$comparisons) {
    r <- run_comparison(Lc, meta, cm, comp_cat, cfg)
    r$table <- cbind(name = names_[rownames(r$table)], r$table)
    res[[cm$name]] <- r$table; summ[[cm$name]] <- r$info; pairs_all[[cm$name]] <- r$pairs
  }
  sheets$summary_univariate <- as.data.frame(data.table::rbindlist(summ, fill = TRUE))
  sheets$age_matched_pairs <- do.call(rbind, pairs_all)
  for (s in summ) log_msg(sprintf("  %s: n %d/%d (в модели %d/%d), пар по возрасту %d, q_model < %.2f: %d, находок: %d",
                                  s$comparison, s$n_case, s$n_ref, s$model_n_case, s$model_n_ref, s$matched_pairs, fdr, s$n_q_model, s$n_findings))
  pa <- feat[feat$kept & feat$presence_absence, c("mz", "rt_min", "annotation", paste0("det_", groups))]
  for (cm in cfg$design$comparisons) {
    s <- comparison_sets(meta, cm)
    pa[[paste0("p_fisher_", cm$name)]] <- vapply(rownames(pa), function(f) {
      a <- sum(!is.na(X1[s$case, f])); b <- sum(!is.na(X1[s$ref, f]))
      stats::fisher.test(matrix(c(a, length(s$case) - a, b, length(s$ref) - b), 2, byrow = TRUE))$p.value }, 0)
  }
  sheets$presence_absence <- pa
  # чувствительность к нормализации
  nc <- list()
  if (length(cfg$preprocess$compare_normalizations) >= 2) {
    m2 <- cfg$preprocess$compare_normalizations[1:2]
    for (cm in cfg$design$comparisons) {
      ra <- run_comparison(norm[[m2[1]]][study, reps], meta, cm, comp_cat, cfg)$table
      rb <- run_comparison(norm[[m2[2]]][study, reps], meta, cm, comp_cat, cfg)$table
      sa <- rownames(ra)[ra$q_model < fdr]; sb <- rownames(rb)[rb$q_model < fdr]
      nc[[cm$name]] <- data.frame(comparison = cm$name, q_A = length(sa), q_B = length(sb), overlap = length(intersect(sa, sb)),
                                  findings_A = sum(ra$finding), findings_B = sum(rb$finding),
                                  findings_overlap = length(intersect(rownames(ra)[ra$finding], rownames(rb)[rb$finding])),
                                  rho_log2FC = stats::cor(ra$beta_model, rb$beta_model, method = "spearman"))
      names(nc[[cm$name]])[2:6] <- c(paste0("q<FDR_", m2), "overlap", paste0("findings_", m2))
      plot_norm_concordance(ra, rb, cm$name, fig(paste0("09_norm_", cm$name)), dpi)
    }
    sheets$normalization_sensitivity <- do.call(rbind, nc)
  }

  # ---------------------------------------------------------------- 5. многомерный анализ
  mv <- list(); curves <- list()
  if (isTRUE(cfg$multivariate$run)) {
    ex <- panel_excluded(comp_cat)
    if (!is.null(meta$contam_index)) {
      rho <- apply(Lc, 2, function(v) stats::cor(v, meta[study, "contam_index"], method = "spearman"))
      feat$rho_contam_index <- NA; feat[reps, "rho_contam_index"] <- rho
      ex <- ex | abs(rho) >= cfg$multivariate$max_abs_rho_contam
    }
    pool <- reps[!ex]
    log_msg(sprintf("Пул для моделей: %d соединений (исключено до CV: %d)", length(pool), sum(ex)))
    for (k in seq_along(cfg$design$comparisons)) {
      cm <- cfg$design$comparisons[[k]]; s <- comparison_sets(meta, cm); pr <- pairs_all[[cm$name]]
      sets <- list(all = s, matched = list(case = pr$case, ref = pr$ref))
      if (!is.null(meta$contam_index)) { zp <- contam_match(meta, s$case, s$ref, cfg); sets$contam_matched <- list(case = zp$case, ref = zp$ref) }
      for (lab in names(sets)) {
        need <- if (lab == "all") cfg$multivariate$cv_folds else cfg$multivariate$min_pairs_cv
        nmin <- min(length(sets[[lab]]$case), length(sets[[lab]]$ref))
        if (nmin < need) { mv[[paste(cm$name, lab)]] <- data.frame(comparison = cm$name, set = lab,
                             model = sprintf("не оценивается: %d < %d на группу", nmin, need)); next }
        e <- evaluate_models(Lc, meta, sets[[lab]]$case, sets[[lab]]$ref, pool, lab, cfg$seed + k, cfg, names_)
        mv[[paste(cm$name, lab)]] <- cbind(comparison = cm$name, e$table)
        curves[[cm$name]] <- c(curves[[cm$name]], e$curves)
      }
    }
    sheets$multivariate_CV <- as.data.frame(data.table::rbindlist(mv, fill = TRUE))
    for (i in seq_len(nrow(sheets$multivariate_CV))) with(sheets$multivariate_CV[i, ], if (!is.na(AUC_mean %||% NA))
      log_msg(sprintf("  %s [%s] %s: AUC %.2f ± %.2f, perm p = %.3g", comparison, set, model, AUC_mean, AUC_sd, perm_p)))
  }

  # ---------------------------------------------------------------- 6. рисунки и отчёт
  find <- as.data.frame(data.table::rbindlist(lapply(names(res), function(n) {
    r <- res[[n]][res[[n]]$finding, ]; if (nrow(r)) cbind(comparison = n, feature = rownames(r), r) }), fill = TRUE))
  sheets$findings <- find
  pngs <- list()
  pngs$qc <- plot_qc_overview(ft$X, meta, feat$RSD_QC_raw[feat$pass_detection], facs[study, , drop = FALSE], fig("01_qc_overview"), dpi)
  pngs$dist <- plot_sample_distributions(log2(X1[ana, ]), L, meta, fig("02_sample_distributions"), dpi)
  full <- reps[colSums(is.na(X1[study, reps, drop = FALSE])) == 0]
  pc <- pca_scores(L[ana, if (length(full) >= 10) full else reps, drop = FALSE])
  pngs$pca <- plot_pca(pc, meta, fig("03_pca"), dpi)
  if (!is.null(ci)) pngs$contam <- plot_contamination(meta, groups, fig("04_contamination_index"), dpi)
  for (cm in cfg$design$comparisons) {
    r <- res[[cm$name]]; n <- cm$name
    pngs[[paste0("volcano_", n)]] <- plot_volcano(r, n, names_, fdr, fig(paste0("05_volcano_", n)), dpi)
    pngs[[paste0("heat_", n)]] <- plot_heatmap(Lc, meta, r, cm$case, cm$ref, names_, cfg$report$heatmap_top, fig(paste0("06_heatmap_", n)))
    top <- rownames(r)[r$finding][order(r$p_model[r$finding])]
    if (!length(top)) top <- rownames(r)[!(r$compound_category %in% EXPOSURE)][order(r$p_model[!(r$compound_category %in% EXPOSURE)])][1:4]
    pngs[[paste0("box_", n)]] <- plot_boxplots(L, meta, utils::head(top, cfg$report$top_boxplots), names_, groups,
                                               paste(n, ": ведущие соединения"), fig(paste0("07_boxplots_", n)), dpi)
    pngs[[paste0("roc_", n)]] <- plot_roc(curves[[n]], paste(n, ": ROC (повторная CV)"), fig(paste0("08_roc_", n)), dpi)
  }
  for (n in names(res)) sheets[[paste0("stats_", n)]] <- res[[n]][order(res[[n]]$p_model), ]
  sheets$samples <- meta
  sheets$features <- feat
  if (!is.null(series)) sheets$contaminant_series <- series
  log_msg(sprintf("Время обработки: %.1f мин", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  sheets <- c(list(log = data.frame(log = log_lines())), sheets)
  write_tables(sheets, out)
  data.table::fwrite(data.frame(feature = colnames(L), t(L)), file.path(out, "tables", "processed_log2_matrix.csv"))
  write_session_info(out)
  man <- write_manifest(out)
  if (isTRUE(cfg$report$html)) {
    sec <- list(
      list(title = "Журнал обработки", html = paste0("<pre>", htmlesc(paste(log_lines(), collapse = "\n")), "</pre>")),
      list(title = "Демография", html = paste(html_table(dm$table), html_table(dm$tests))),
      list(title = "Контроль качества", html = paste(img_tag(pngs$qc, "TIC, RSD в QC, факторы разведения"),
                                                    img_tag(pngs$dist, "Распределения интенсивностей по пробам"),
                                                    img_tag(pngs$pca, "PCA"), img_tag(pngs$contam, "Индекс загрязнения"))),
      list(title = "Одномерная статистика", html = paste(html_table(sheets$summary_univariate),
             paste(vapply(names(res), function(n) paste(img_tag(pngs[[paste0("volcano_", n)]], paste("Вулкан:", n)),
                                                        img_tag(pngs[[paste0("heat_", n)]], paste("Тепловая карта:", n)),
                                                        img_tag(pngs[[paste0("box_", n)]], paste("Диаграммы размаха:", n))), ""), collapse = ""))),
      list(title = "Классификация (кросс-валидация)", html = paste(html_table(sheets$multivariate_CV),
             paste(vapply(names(res), function(n) img_tag(pngs[[paste0("roc_", n)]], paste("ROC:", n)), ""), collapse = ""))),
      list(title = "Находки", html = html_table(utils::head(find[, intersect(c("comparison", "name", "beta_model", "q_model", "p_matched", "AUC", "compound_category"), names(find %||% data.frame()))], 100))),
      list(title = "Контрольные суммы таблиц (md5)", html = html_table(man)))
    write_html_report(out, cfg, sec)
  }
  log_msg("Готово → ", normalizePath(out))
  invisible(list(meta = meta, L = L, Lc = Lc, res = res, features = feat, sheets = sheets))
}

if (sys.nframe() == 0L) main()

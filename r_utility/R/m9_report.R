# Модуль 9. База метаболомных признаков (CSV с дозаписью) и отчёт HTML (R Markdown; без pandoc —
# встроенный генератор HTML) со всеми таблицами, рисунками, журналом и sessionInfo().

db_append <- function(df, path) {
  if (file.exists(path)) {
    h <- names(data.table::fread(path, nrows = 0))
    if (!identical(h, names(df))) stop("Структура базы ", path, " не совпадает с текущей версией — укажите другую папку базы")
    data.table::fwrite(df, path, append = TRUE)
  } else data.table::fwrite(df, path, bom = TRUE)
}

module9_report <- function(s) {
  out <- module_dir(s, 9); d <- load_module(s, 7)
  run_id <- .ue$run_id %||% format(Sys.time(), "%Y%m%d_%H%M%S")
  # ---------------- база признаков
  db <- project_path(s, s$report$database); dir.create(db, showWarnings = FALSE, recursive = TRUE)
  m3 <- load_module(s, 3); m4 <- load_module(s, 4)
  u <- d$sel$uvf
  prof <- s$report$database_profile
  ids <- intersect(if (prof == "strict") m3$strict_ids %||% colnames(d$L) else m3$standard_ids %||% colnames(d$L), colnames(d$L))
  if (!length(ids)) {                           # например, строгий профиль пуст при большой доле пропусков
    log_msg("  ВНИМАНИЕ: в профиле базы признаков нет признаков — в базу записывается профиль анализа")
    prof <- paste0(prof, "→", s$filter$analysis_profile); ids <- colnames(d$L)
  }
  ft <- data.frame(run_id = run_id, date = format(Sys.Date()), project = s$project$name, db_profile = prof, id = ids,
                   d$info[match(ids, d$info$id), c("mz", "rt_min", "name")],
                   rsd_qc = m3$filter$rsd_qc[match(ids, m3$filter$id)] %||% NA,
                   p_uvf = u$p[match(ids, u$id)], p_adj = u$p_adj[match(ids, u$id)],
                   log2FC = if (!is.null(u$log2FC)) u$log2FC[match(ids, u$id)] else NA,
                   selected_uvf = ids %in% d$sel$f_uvf,
                   n_models = if (!is.null(d$sel$freq)) d$sel$freq$n_models[match(ids, d$sel$freq$id)] else NA,
                   selected_spf = ids %in% d$sel$f_spf, selected_rfe = ids %in% d$sel$f_rfe,
                   nested_share = if (!is.null(d$nested)) d$nested$frequency$share_of_folds[match(ids, d$nested$frequency$id)] else NA)
  X4 <- m4$X[, ids, drop = FALSE]
  inten <- data.frame(run_id = run_id, sample = rep(rownames(X4), ncol(X4)), group = rep(m4$meta$group, ncol(X4)),
                      id = rep(colnames(X4), each = nrow(X4)),
                      mz = rep(d$info$mz[match(colnames(X4), d$info$id)], each = nrow(X4)),
                      rt_min = rep(d$info$rt_min[match(colnames(X4), d$info$id)], each = nrow(X4)),
                      intensity_raw = as.vector(X4), log2_corrected = as.vector(d$L[rownames(X4), colnames(X4)]))
  runs <- data.frame(run_id = run_id, date = format(Sys.time(), "%Y-%m-%d %H:%M"), project = s$project$name, version = UROMETABR_VERSION,
                     input = basename(if (s$input$type == "table") s$input$table_file else s$input$raw_dir),
                     n_samples = nrow(d$L), groups = paste(sprintf("%s:%d", names(table(d$meta$group)), table(d$meta$group)), collapse = " "),
                     n_features = ncol(d$L), db_profile = prof, n_db_features = length(ids), n_uvf = length(d$sel$f_uvf), n_spf = length(d$sel$f_spf), n_rfe = length(d$sel$f_rfe),
                     p_adjust = s$select$p_adjust, seed = s$project$seed)
  db_append(runs, file.path(db, "runs_db.csv")); db_append(ft, file.path(db, "features_db.csv")); db_append(inten, file.path(db, "intensities_db.csv"))
  log_msg(sprintf("  База признаков (%s, %s профиль): добавлено %d признаков, %d значений интенсивности (прогон %s)", db,
                  if (prof == "strict") "строгий" else "стандартный", nrow(ft), nrow(inten), run_id))
  # ---------------- отчёт
  sec <- report_sections(s, run_id)
  saveRDS(sec, file.path(out, "report_sections.rds"))
  f <- file.path(out, "report.html")
  rendered <- FALSE
  if (isTRUE(s$report$rmarkdown) && requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    tpl <- file.path(getOption("urometabr.root", "."), "templates", "report.Rmd")
    rendered <- tryCatch({
      rmarkdown::render(tpl, output_file = "report.html", output_dir = out, intermediates_dir = tempdir(),
                        params = list(sections = file.path(out, "report_sections.rds")), quiet = TRUE, envir = new.env())
      TRUE }, error = function(e) { log_msg("  R Markdown недоступен (", conditionMessage(e), ") — встроенный генератор HTML"); FALSE })
  }
  if (!rendered) write_html(sec, f)
  log_msg("  Отчёт: ", f, if (rendered) " (R Markdown)" else " (встроенный генератор)")
  save_module(s, 9, list(report = f, database = db, run_id = run_id))
}

# Содержание отчёта: список разделов (текст, таблицы, рисунки) — общий для R Markdown и встроенного HTML
report_sections <- function(s, run_id) {
  P <- function(id, f) file.path(s$project$dir, MODULES$dir[MODULES$id == id], f)
  rd <- function(id, f) if (file.exists(P(id, f))) as.data.frame(data.table::fread(P(id, f), encoding = "UTF-8")) else NULL
  img <- function(id, f) { p <- P(id, paste0(f, ".png")); if (file.exists(p)) p else NULL }
  logs <- file.path(s$project$dir, "logs")
  lastf <- function(pat) { f <- sort(list.files(logs, pattern = pat, full.names = TRUE)); if (length(f)) utils::tail(f, 1) else NULL }
  m7 <- load_module(s, 7)
  sec <- list(
    list(title = "Сведения о прогоне", text = c(sprintf("%s. Версия %s (%s).", PROGRAM_NAME, UROMETABR_VERSION, UROMETABR_DATE),
           sprintf("Проект: %s; папка: %s; прогон %s; seed %s.", s$project$name, s$project$dir, run_id, s$project$seed)),
         tables = list("Время и память по модулям" = if (!is.null(lastf("^modules_"))) as.data.frame(data.table::fread(lastf("^modules_"), encoding = "UTF-8")))),
    list(title = "0. Входные данные", tables = list("Проверки" = rd(0, "input_check.csv"))),
    list(title = "2. Контроль качества", tables = list("Сводка" = rd(2, "qc_summary.csv"), "Пропуски по группам" = rd(2, "missing_by_group.csv")),
         images = c(img(2, "pca_qc"), img(2, "tic_by_order"), img(2, "rsd_qc"))),
    list(title = "3. Фильтр устойчивых признаков", text = "Стандартный профиль — для анализа (модули 4–8); строгий — для базы признаков (модуль 9).",
         tables = list("Шаги фильтрации по профилям" = rd(3, "filter_steps.csv"))),
    list(title = "5. Заполнение пропусков", tables = list(rd(5, "impute_info.csv"))),
    list(title = "6. Коррекция сигнала", tables = list(rd(6, "correct_info.csv")), images = img(6, "distribution_by_sample")),
    list(title = "7. Статистика и отбор признаков",
         text = "Средняя AUROC по этапам — описательная величина (признаки отобраны по тем же данным). Несмещённая оценка — вложенная кросс-валидация.",
         tables = list("Число признаков и средняя AUROC по этапам" = rd(7, "auroc_steps.csv"), "Точность моделей (кросс-валидация)" = rd(7, "model_accuracy.csv"),
                       "Вложенная кросс-валидация" = rd(7, "nested_cv.csv"),
                       "Признаки после РОП" = { x <- rd(7, "uvf_stats.csv"); if (!is.null(x)) x[x$id %in% m7$sel$f_rfe, ] },
                       "ОСФ: 30 признаков с наименьшими p" = utils::head(rd(7, "uvf_stats.csv"), 30),
                       "Модель с поправкой на пол и возраст (30 первых)" = utils::head(rd(7, "covariate_model.csv"), 30))),
    list(title = "8. Визуализация", images = c(img(8, "pca_stages"), img(8, "hca_stages"), img(8, "heatmap"), img(8, "boxplots_metabolites"), img(8, "plsda_vip"))),
    list(title = "9. База признаков", text = sprintf("Записи добавлены в %s (runs_db.csv, features_db.csv, intensities_db.csv).", project_path(s, s$report$database))),
    list(title = "Журнал прогона", pre = if (!is.null(.ue$log_file)) readLines(.ue$log_file, encoding = "UTF-8")),
    list(title = "sessionInfo()", pre = utils::capture.output(utils::sessionInfo())))
  sec
}

html_esc <- function(x) { x <- gsub("&", "&amp;", x); x <- gsub("<", "&lt;", x); gsub(">", "&gt;", x) }
html_tab <- function(d) {
  if (is.null(d) || !nrow(d)) return("")
  f <- function(v) if (is.numeric(v)) formatC(v, digits = 4, format = "g") else html_esc(as.character(v))
  rows <- apply(as.data.frame(lapply(d, f)), 1, function(r) paste0("<tr><td>", paste(r, collapse = "</td><td>"), "</td></tr>"))
  paste0("<table><tr><th>", paste(html_esc(names(d)), collapse = "</th><th>"), "</th></tr>", paste(rows, collapse = ""), "</table>")
}
write_html <- function(sec, f) {
  body <- vapply(sec, function(x) {
    tb <- if (length(x$tables)) paste(vapply(seq_along(x$tables), function(i) {
      nm <- names(x$tables)[i]; paste0(if (!is.null(nm) && nzchar(nm)) paste0("<h4>", html_esc(nm), "</h4>") else "", html_tab(x$tables[[i]])) }, ""), collapse = "") else ""
    im <- paste(vapply(x$images %||% character(0), function(p)
      sprintf("<img src=\"data:image/png;base64,%s\"/>", jsonlite::base64_enc(readBin(p, "raw", file.info(p)$size))), ""), collapse = "")
    paste0("<h2>", html_esc(x$title), "</h2>", paste0("<p>", html_esc(x$text %||% character(0)), "</p>", collapse = ""), tb, im,
           if (length(x$pre)) paste0("<pre>", html_esc(paste(x$pre, collapse = "\n")), "</pre>") else "")
  }, "")
  css <- "body{font-family:Arial,sans-serif;max-width:1100px;margin:20px auto;padding:0 16px}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccc;padding:2px 6px}th{background:#eee}img{max-width:100%;margin:8px 0}pre{background:#f6f6f6;padding:8px;font-size:11px;white-space:pre-wrap}"
  writeLines(c("<!doctype html><html lang=\"ru\"><head><meta charset=\"utf-8\"><title>Отчёт UroMetabR</title><style>", css, "</style></head><body>",
               paste0("<h1>Отчёт о прогоне UroMetabR ", UROMETABR_VERSION, "</h1>"), body, "</body></html>"), f, useBytes = TRUE)
}

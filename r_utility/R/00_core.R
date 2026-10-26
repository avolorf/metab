# Ядро утилиты: версия, настройки, папка проекта, журнал, запуск модулей с замером времени и памяти.

UROMETABR_VERSION <- "1.0.0"
UROMETABR_DATE <- "2026-10-26"
PROGRAM_NAME <- paste("Программный инструмент для неинвазивной вспомогательной диагностики рака",
                      "мочеполовой системы на основе анализа мочи методом ВЭЖХ-МС")

MODULES <- data.frame(
  id = 0:9,
  dir = c("00_input", "01_peaktable", "02_qc", "03_filter", "04_samples", "05_impute", "06_correct",
          "07_select", "08_visual", "09_report"),
  title = c("Загрузка и проверка входных данных", "Оптимизация параметров и таблица пиков (xcms)",
            "Контроль качества", "Фильтр устойчивых признаков", "Удаление QC, усреднение повторов, сортировка",
            "Заполнение пропусков", "Коррекция сигнала и логарифмирование", "Статистика и отбор признаков",
            "Визуализация", "База признаков и отчёт"),
  stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ------------------------------------------------------------------ настройки
default_settings <- function() list(
  project = list(name = "Проект", seed = 20261026, cores = 1),
  input = list(type = "table",                       # table (CSV/TSV/TXT, экспорт MS-DIAL) | raw (mzML/mzXML/CDF)
               table_file = "", raw_dir = "", metadata_file = "",
               table_format = "auto",                # auto | wide | msdial
               id_column = ""),                      # пусто — первый столбец
  xcms = list(ipo = FALSE, ppm = 5, min_peakwidth = 5, max_peakwidth = 30, snthresh = 10, noise = 1000,
              prefilter_k = 3, prefilter_i = 1000, bw = 5, min_fraction = 0.5, mzwid = 0.015,
              obiwarp_binsize = 0.5),
  filter = list(min_qc_detect = 0.8, min_group_detect = 0.5, max_qc_rsd = 0.30, blank_ratio = 5,
                strict = list(min_all_detect = 1.0, min_qc_detect = 0.8, max_qc_rsd = 0.15),   # профиль базы признаков
                analysis_profile = "standard",                                                  # standard | strict
                target_min = 300, target_max = 800),
  samples = list(average_replicates = TRUE, groups = list()),   # пусто — все группы, кроме QC и бланков
  impute = list(method = "rf", maxiter = 10, ntree = 100),      # rf | missForest | halfmin
  correct = list(normalization = "none", log2 = TRUE, method = "eigenms", eigen_perm = 200, eigen_alpha = 0.05),
  select = list(p_adjust = "holm",                              # holm | BH | none | stage1
                alpha = 0.05, models = c("knn", "svmRadial", "pls", "rf"), cv_folds = 5, cv_repeats = 3,
                tune_length = 5, top_n = 50, min_freq = 2, rfe_method = "nb", rfe_folds = 5, rfe_repeats = 3,
                compare_stages = TRUE, nested_cv = TRUE, nested_folds = 5, nested_repeats = 2,
                nested_input = "imputed"),                       # imputed (без EigenMS) | corrected
  visual = list(dpi = 600, formats = c("png", "jpeg"), heatmap_top = 50, boxplot_top = 12, vip_top = 20,
                hca_distance = "manhattan", hca_method = "ward.D2"),
  report = list(database = "features_db", database_profile = "strict", rmarkdown = TRUE),   # strict | standard
  run = list(from = 0, to = 9))

merge_settings <- function(base, over) {
  for (k in names(over)) {
    if (is.list(base[[k]]) && is.list(over[[k]]) && length(over[[k]]) && !is.null(names(over[[k]])))
      base[[k]] <- merge_settings(base[[k]], over[[k]])
    else base[k] <- list(over[[k]])
  }
  base
}

read_settings <- function(path) {
  if (!file.exists(path)) stop("Файл настроек не найден: ", path)
  s <- merge_settings(default_settings(), yaml::read_yaml(path, fileEncoding = "UTF-8") %||% list())
  s$project$dir <- normalizePath(dirname(path), winslash = "/")
  s$settings_file <- normalizePath(path, winslash = "/")
  check_settings(s)
  s
}

write_settings <- function(s, path) {
  s$project$dir <- NULL; s$settings_file <- NULL
  s <- rapply(s, function(x) if (all(x == round(x)) && all(abs(x) < 2^31)) as.integer(x) else x, classes = "numeric", how = "replace")
  yaml::write_yaml(s, path, fileEncoding = "UTF-8")
  invisible(path)
}

# Путь из настроек: абсолютный — как есть, относительный — от папки проекта
project_path <- function(s, p) {
  if (is.null(p) || !nzchar(p)) return("")
  if (grepl("^(/|[A-Za-z]:|\\\\\\\\)", p)) p else file.path(s$project$dir, p)
}

check_settings <- function(s) {
  e <- character(0)
  if (!s$input$type %in% c("table", "raw")) e <- c(e, "input$type: table или raw")
  if (!s$impute$method %in% c("rf", "missForest", "halfmin")) e <- c(e, "impute$method: rf, missForest или halfmin")
  if (!s$correct$method %in% c("eigenms", "none")) e <- c(e, "correct$method: eigenms или none")
  if (!s$correct$normalization %in% c("none", "PQN")) e <- c(e, "correct$normalization: none или PQN")
  if (!s$select$p_adjust %in% c("holm", "BH", "none", "stage1")) e <- c(e, "select$p_adjust: holm, BH, none или stage1")
  if (!s$select$rfe_method %in% c("nb", "rf")) e <- c(e, "select$rfe_method: nb или rf")
  if (!s$select$nested_input %in% c("imputed", "corrected")) e <- c(e, "select$nested_input: imputed или corrected")
  if (!s$filter$analysis_profile %in% c("standard", "strict")) e <- c(e, "filter$analysis_profile: standard или strict")
  if (!s$report$database_profile %in% c("standard", "strict")) e <- c(e, "report$database_profile: standard или strict")
  if (s$filter$max_qc_rsd <= 0 || s$filter$min_group_detect < 0 || s$filter$min_group_detect > 1) e <- c(e, "пороги фильтра вне допустимых значений")
  if (length(e)) stop("Ошибки в настройках:\n  ", paste(e, collapse = "\n  "))
  invisible(TRUE)
}

# ------------------------------------------------------------------ журнал
.ue <- new.env()
.ue$log_file <- NULL

log_start <- function(s) {
  d <- file.path(s$project$dir, "logs"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  .ue$run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  .ue$log_file <- file.path(d, paste0("run_", .ue$run_id, ".log"))
  file.create(.ue$log_file)
  .ue$run_id
}

log_msg <- function(...) {
  s <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
  cat(s, "\n", sep = "")
  if (!is.null(.ue$log_file)) cat(s, "\n", sep = "", file = .ue$log_file, append = TRUE)
  invisible(s)
}

# ------------------------------------------------------------------ результаты модулей
module_dir <- function(s, id) {
  d <- file.path(s$project$dir, MODULES$dir[MODULES$id == id])
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}
save_module <- function(s, id, obj) saveRDS(obj, file.path(module_dir(s, id), "result.rds"))
load_module <- function(s, id) {
  f <- file.path(s$project$dir, MODULES$dir[MODULES$id == id], "result.rds")
  if (!file.exists(f)) stop(sprintf("Нет результата модуля %d (%s): запустите его раньше", id, MODULES$title[MODULES$id == id]))
  readRDS(f)
}
# Последний доступный результат с номером < id (для таблиц модуль 1 пропускается)
load_previous <- function(s, id) {
  for (k in rev(seq_len(id) - 1)) {
    f <- file.path(s$project$dir, MODULES$dir[MODULES$id == k], "result.rds")
    if (file.exists(f)) return(readRDS(f))
  }
  stop("Нет результатов предыдущих модулей — запустите модуль 0")
}
write_csv <- function(x, path) data.table::fwrite(x, path, bom = TRUE)   # BOM — корректная кириллица в Excel

# ------------------------------------------------------------------ запуск конвейера
MODULE_FUN <- c("module0_input", "module1_peaktable", "module2_qc", "module3_filter", "module4_samples",
                "module5_impute", "module6_correct", "module7_select", "module8_visual", "module9_report")

run_pipeline <- function(s, from = s$run$from, to = s$run$to) {
  if (is.character(s)) s <- read_settings(s)
  log_start(s)
  log_msg("UroMetabR ", UROMETABR_VERSION, " (", UROMETABR_DATE, "); R ", R.version$major, ".", R.version$minor,
          "; ", utils::sessionInfo()$running)
  log_msg("Проект: ", s$project$dir, "; модули ", from, "–", to, "; seed ", s$project$seed)
  file.copy(s$settings_file %||% "", file.path(s$project$dir, "logs", paste0("settings_", .ue$run_id, ".yaml")))
  status <- data.frame()
  for (id in from:to) {
    if (id == 1 && s$input$type == "table") {
      log_msg("Модуль 1 пропущен: на входе готовая таблица признаков")
      status <- rbind(status, data.frame(module = 1, title = MODULES$title[2], status = "пропущен", time_s = 0, peak_ram_mb = NA))
      next
    }
    log_msg(sprintf("=== Модуль %d. %s", id, MODULES$title[MODULES$id == id]))
    set.seed(s$project$seed + id)                   # единый seed: результат модуля не зависит от предыдущих запусков
    invisible(gc(reset = TRUE)); t0 <- proc.time()[["elapsed"]]
    ok <- tryCatch({ get(MODULE_FUN[id + 1])(s); TRUE }, error = function(e) { log_msg("ОШИБКА: ", conditionMessage(e)); FALSE })
    tm <- proc.time()[["elapsed"]] - t0
    mem <- sum(gc()[, 6])
    log_msg(sprintf("Модуль %d: %s, %.1f с, пиковая память R %.0f МБ", id, if (ok) "выполнен" else "ошибка", tm, mem))
    status <- rbind(status, data.frame(module = id, title = MODULES$title[id + 1], status = if (ok) "выполнен" else "ошибка",
                                       time_s = round(tm, 1), peak_ram_mb = round(mem)))
    write_csv(status, file.path(s$project$dir, "logs", paste0("modules_", .ue$run_id, ".csv")))
    if (!ok) { write_session_info(s); log_msg("Обработка остановлена на модуле ", id); return(invisible(status)) }
  }
  write_session_info(s)
  log_msg("Готово. Отчёт: ", file.path(s$project$dir, "09_report", "report.html"))
  cat("UROMETABR_DONE\n", file = .ue$log_file, append = TRUE)
  invisible(status)
}

write_session_info <- function(s) {
  f <- file.path(s$project$dir, "logs", paste0("session_info_", .ue$run_id, ".txt"))
  writeLines(c(paste("UroMetabR", UROMETABR_VERSION, UROMETABR_DATE), utils::capture.output(utils::sessionInfo())), f)
  f
}

# ------------------------------------------------------------------ параллельные вычисления
# Linux/macOS — fork; Windows — PSOCK-кластер с подгрузкой модулей утилиты
par_lapply <- function(X, FUN, cores = 1) {
  if (cores <= 1 || length(X) < 2) return(lapply(X, FUN))
  if (.Platform$OS.type == "unix") return(parallel::mclapply(X, FUN, mc.cores = cores))
  cl <- parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl))
  files <- list.files(getOption("urometabr.dir"), pattern = "\\.R$", full.names = TRUE)
  parallel::clusterCall(cl, function(f) for (x in f) source(x, encoding = "UTF-8"), files)
  parallel::parLapply(cl, X, FUN)
}

# Зерна для caret: воспроизводимость и в параллельном режиме
caret_seeds <- function(n_resamples, n_tune, seed) {
  set.seed(seed)
  c(lapply(seq_len(n_resamples), function(i) sample.int(1e6, n_tune)), list(sample.int(1e6, 1)))
}

register_cores <- function(cores) {
  if (cores > 1 && requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(caret)))   # модели caret ищут функции в caret
    doParallel::registerDoParallel(cl); return(cl)
  }
  foreach::registerDoSEQ(); NULL
}
stop_cores <- function(cl) { if (!is.null(cl)) parallel::stopCluster(cl); foreach::registerDoSEQ() }

# ------------------------------------------------------------------ графика
GROUP_PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
group_colors <- function(g) {
  u <- sort(unique(as.character(g))); u <- c(setdiff(u, "QC"), intersect(u, "QC"))
  col <- stats::setNames(GROUP_PAL[seq_along(u)], u); if ("QC" %in% u) col["QC"] <- "#7F7F7F"
  col
}
save_fig <- function(s, p, dir, name, w = 17, h = 12) {
  out <- character(0)
  for (fmt in s$visual$formats) {
    f <- file.path(dir, paste0(name, ".", if (fmt == "jpeg") "jpg" else fmt))
    ggplot2::ggsave(f, p, width = w, height = h, units = "cm", dpi = s$visual$dpi, bg = "white")
    out <- c(out, f)
  }
  invisible(out)
}
# Графики не-ggplot (дендрограммы, тепловые карты): одна функция рисования — несколько форматов
save_base_fig <- function(s, draw, dir, name, w = 17, h = 12) {
  for (fmt in s$visual$formats) {
    f <- file.path(dir, paste0(name, ".", if (fmt == "jpeg") "jpg" else fmt))
    dev <- if (fmt == "jpeg") grDevices::jpeg else grDevices::png
    args <- list(f, width = w, height = h, units = "cm", res = s$visual$dpi)
    if (capabilities("cairo")) args$type <- "cairo"
    if (fmt == "jpeg") args$quality <- 95
    do.call(dev, args); try(draw()); grDevices::dev.off()
  }
}
theme_um <- function() ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top")

# UroMetabR — общие функции: журнал, параметры, проверки, параллельные вычисления.

UROMETABR_VERSION <- "1.0.0"

.log_env <- new.env()
.log_env$lines <- character(0)

log_msg <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  .log_env$lines <- c(.log_env$lines, s)
  invisible(s)
}

log_reset <- function() .log_env$lines <- character(0)
log_lines <- function() .log_env$lines

`%||%` <- function(a, b) if (is.null(a)) b else a

# Параметры по умолчанию; файл config.yaml переопределяет только указанные поля.
default_config <- function() list(
  input = list(data = NULL, format = "wide", metadata = NULL, id_column = "Sample",
               meta_columns = c("Sample", "Group", "Diag"),
               contaminant_series = NULL, exposure_keywords = NULL, mzml_dir = NULL),
  output = "results_R",
  seed = 20261013,
  cores = 1,
  design = list(groups = c("Control", "BC", "PC"),
                qc_pattern = "^QC", blank_pattern = "^(blank|бланк|blk)", qc_every = 20,
                comparisons = list(list(name = "BC_vs_Control", case = "BC", ref = "Control", sex = NULL)),
                covariates = c("sex", "age", "contam_index"),
                match_caliper_years = 5, match_within_sex = TRUE),
  preprocess = list(min_qc_detect = 0.8, min_group_detect = 0.5, max_qc_rsd = 0.30, blank_ratio = 5,
                    drift = "none", normalization = "PQN", compare_normalizations = c("PQN", "MSTUS"),
                    low_dilution_factor = 0.2, presence_absence = c(0.5, 0.8)),
  grouping = list(rt_tol = 0.02, min_corr = 0.8, ppm = 5),
  stats = list(fdr = 0.05, confirm_p = 0.05, ctrl_index_p = 0.05),
  multivariate = list(cv_folds = 5, cv_repeats = 20, panel_k = 5, n_perm = 1000,
                      max_abs_rho_contam = 0.5, min_pairs_cv = 10, run = TRUE),
  report = list(fig_dpi = 300, top_boxplots = 12, heatmap_top = 50, html = TRUE)
)

merge_list <- function(base, over) {
  for (k in names(over)) {
    if (is.list(base[[k]]) && is.list(over[[k]]) && !is.null(names(over[[k]])) && k != "comparisons")
      base[[k]] <- merge_list(base[[k]], over[[k]])
    else base[k] <- list(over[[k]])
  }
  base
}

read_config <- function(path = NULL, quick = FALSE) {
  cfg <- default_config()
  if (!is.null(path)) {
    if (!file.exists(path)) stop("Файл параметров не найден: ", path)
    cfg <- merge_list(cfg, yaml::read_yaml(path))
    base <- dirname(normalizePath(path))
    # относительные пути в config.yaml считаются от каталога, где лежит config.yaml
    for (k in c("data", "metadata", "contaminant_series", "exposure_keywords", "mzml_dir")) {
      p <- cfg$input[[k]]
      if (!is.null(p) && nzchar(p) && !grepl("^(/|[A-Za-z]:)", p)) cfg$input[[k]] <- file.path(base, p)
    }
    if (!grepl("^(/|[A-Za-z]:)", cfg$output)) cfg$output <- file.path(base, cfg$output)
  }
  if (quick) { cfg$multivariate$n_perm <- 20; cfg$multivariate$cv_repeats <- 3 }
  validate_config(cfg)
  cfg
}

validate_config <- function(cfg) {
  p <- cfg$preprocess
  stopifnot(p$min_qc_detect >= 0, p$min_qc_detect <= 1, p$min_group_detect >= 0, p$min_group_detect <= 1,
            p$max_qc_rsd > 0, p$blank_ratio > 0)
  if (!p$drift %in% c("none", "lowess", "auto")) stop("preprocess$drift: none | lowess | auto")
  if (!p$normalization %in% c("PQN", "MSTUS", "none")) stop("preprocess$normalization: PQN | MSTUS | none")
  if (!cfg$input$format %in% c("wide", "msdial", "xcms")) stop("input$format: wide | msdial | xcms")
  for (cm in cfg$design$comparisons)
    if (!all(c(cm$case, cm$ref) %in% cfg$design$groups)) stop("Сравнение ", cm$name, ": группы вне design$groups")
  invisible(TRUE)
}

# Параллельное применение функции; на Windows — PSOCK-кластер с подгрузкой модулей утилиты.
par_lapply <- function(X, FUN, cores = 1, ...) {
  if (cores <= 1 || length(X) < 2) return(lapply(X, FUN, ...))
  if (.Platform$OS.type == "unix") return(parallel::mclapply(X, FUN, ..., mc.cores = cores))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  files <- list.files(module_dir(), pattern = "\\.R$", full.names = TRUE)
  parallel::clusterCall(cl, function(f) for (x in f) source(x, encoding = "UTF-8"), files)
  parallel::parLapply(cl, X, FUN, ...)
}

module_dir <- function() getOption("urometabr.dir", "R")

# Контрольная сумма таблиц результатов (проверка воспроизводимости)
md5_files <- function(files) vapply(files, function(f) digest::digest(file = f, algo = "md5"), "")

col_medians <- function(M) apply(M, 2, stats::median, na.rm = TRUE)
row_medians <- function(M) apply(M, 1, stats::median, na.rm = TRUE)

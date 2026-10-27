#!/usr/bin/env Rscript
# Проверка модуля 1 (xcms) на открытых наборах Bioconductor: faahKO (12 файлов netCDF: 6 KO, 6 WT) и msdata
# (файлы mzML прибора Sciex). Одна команда из папки утилиты:
#   Rscript tests/module1_check.R --out=module1_check
# Устанавливает при необходимости BiocManager, xcms, MSnbase, faahKO, msdata (нужен интернет), конвертирует
# faahKO в mzML и mzXML (MSnbase::writeMSData) и прогоняет утилиту на трёх форматах. Итог — module1_check/summary.csv,
# журналы и sessionInfo() — в module1_check/<проект>/logs.
if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows")
  for (loc in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
extra <- Sys.getenv("UROMETABR_EXTRA_LIB")              # необязательно: дополнительная библиотека пакетов (в конце .libPaths)
if (nzchar(extra)) .libPaths(c(.libPaths(), strsplit(extra, .Platform$path.sep)[[1]]))
args <- commandArgs(trailingOnly = TRUE)
OUT <- normalizePath(sub("^--out=", "", c(grep("^--out=", args, value = TRUE), "--out=module1_check")[1]), mustWork = FALSE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(normalizePath(self)), ".."))
repos <- c(CRAN = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
need <- c("xcms", "MSnbase", "faahKO", "msdata")
miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE)
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8")
if (length(missing_packages())) install_missing()
load_urometabr(root)
source(file.path(root, "tests", "make_test_data.R"), encoding = "UTF-8")      # new_project()
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- data.frame(check = character(0), result = character(0), detail = character(0))
add <- function(check, ok, detail) { S[nrow(S) + 1, ] <<- list(check, if (isTRUE(ok)) "пройден" else "не пройден", detail)
  data.table::fwrite(S, file.path(OUT, "summary.csv"), bom = TRUE); cat(sprintf("[%s] %s: %s\n", S$result[nrow(S)], check, detail)) }
add("Версии", TRUE, sprintf("R %s.%s; xcms %s; MSnbase %s", R.version$major, R.version$minor, packageVersion("xcms"), packageVersion("MSnbase")))

# ---- faahKO: исходные CDF и конвертация в mzML / mzXML
cdf <- list.files(system.file("cdf", package = "faahKO"), recursive = TRUE, full.names = TRUE, pattern = "\\.CDF$")
md <- utils::read.csv(file.path(root, "tests", "module1", "metadata_faahKO.csv"))
dirs <- c(CDF = file.path(OUT, "faahKO_cdf"), mzML = file.path(OUT, "faahKO_mzML"), mzXML = file.path(OUT, "faahKO_mzXML"))
for (d in dirs) { unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE) }   # повторный запуск: файлы конвертации заново
file.copy(cdf, dirs["CDF"], overwrite = TRUE)
raw <- MSnbase::readMSData(cdf, mode = "onDisk")
nm <- tools::file_path_sans_ext(basename(cdf))
MSnbase::writeMSData(raw, file.path(dirs["mzML"], paste0(nm, ".mzML")), outformat = "mzml", copy = FALSE)
MSnbase::writeMSData(raw, file.path(dirs["mzXML"], paste0(nm, ".mzXML")), outformat = "mzxml", copy = FALSE)
add("faahKO: файлы", length(cdf) == 12, sprintf("%d CDF → mzML и mzXML", length(cdf)))

nfeat <- list()
for (fmt in names(dirs)) {
  p <- file.path(OUT, paste0("project_faahKO_", fmt)); unlink(p, recursive = TRUE); dir.create(p)
  file.copy(file.path(root, "tests", "module1", "metadata_faahKO.csv"), file.path(p, "metadata.csv"))
  f <- new_project(p, "", "metadata.csv", type = "raw", raw_dir = normalizePath(dirs[fmt]),
                   xcms = list(ppm = 25, min_peakwidth = 20, max_peakwidth = 80, snthresh = 10, noise = 1000, prefilter_i = 5000, bw = 30, mzwid = 0.1),
                   select = list(cv_folds = 3, rfe_folds = 3, nested_cv = FALSE, compare_stages = FALSE), visual = list(dpi = 150))
  to <- if (fmt == "mzML") 9 else 3
  t0 <- proc.time()[["elapsed"]]; st <- run_pipeline(read_settings(f), 0, to); tm <- proc.time()[["elapsed"]] - t0
  ok <- all(st$status == "выполнен")
  pt <- file.path(p, "01_peaktable", "peak_table.csv")
  n <- if (file.exists(pt)) ncol(data.table::fread(pt, nrows = 1)) - 1 else NA; nfeat[[fmt]] <- n
  add(sprintf("faahKO %s: модули 0–%d", fmt, to), ok,
      sprintf("признаков в таблице пиков %s; время %.1f мин; модуль 1 — %s с, пиковая память R %s МБ", n, tm / 60,
              st$time_s[st$module == 1], st$peak_ram_mb[st$module == 1]))
}
add("faahKO: одинаковое число признаков для CDF, mzML, mzXML", !anyNA(unlist(nfeat)) && length(unique(unlist(nfeat))) == 1, paste(names(nfeat), unlist(nfeat), collapse = "; "))

# ---- msdata: проверка чтения mzML (Sciex) модулем 0 и MSnbase
ms <- list.files(system.file(package = "msdata"), recursive = TRUE, full.names = TRUE, pattern = "\\.(mzML|mzXML)$")
sc <- ms[grepl("sciex", ms, ignore.case = TRUE)]
for (f in c(sc, ms[grepl("mzXML$", ms)][1])) if (!is.na(f)) {
  chk <- check_raw_file(f); r <- tryCatch(MSnbase::readMSData(f, mode = "onDisk"), error = function(e) NULL)
  add(paste("msdata:", basename(f)), chk$ok && !is.null(r), sprintf("модуль 0: %s; MSnbase: %s спектров", chk$msg, if (is.null(r)) "ошибка" else length(r)))
}
writeLines(utils::capture.output(utils::sessionInfo()), file.path(OUT, "session_info_module1.txt"))
cat("\nИтог:", file.path(OUT, "summary.csv"), "\n")

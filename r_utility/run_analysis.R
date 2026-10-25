#!/usr/bin/env Rscript
# Запуск без графического интерфейса:
#   Rscript run_analysis.R путь/к/settings.yaml [--from=N] [--to=M]
# settings.yaml лежит в папке проекта; результаты модулей 0–9 сохраняются в подпапки 00_input … 09_report.

if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows")      # кириллица в модулях: UTF-8
  for (loc in c("C.UTF-8", "en_US.UTF-8", "ru_RU.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
args <- commandArgs(trailingOnly = TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- if (!is.na(self)) dirname(normalizePath(self)) else getwd()
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8")
load_urometabr(root)
if (length(missing_packages())) install_missing()
cfg <- args[!startsWith(args, "--")][1]
if (is.na(cfg)) stop("Укажите файл настроек: Rscript run_analysis.R папка_проекта/settings.yaml")
s <- read_settings(cfg)
opt <- function(k, d) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) as.integer(sub(".*=", "", a)) else d }
st <- run_pipeline(s, opt("from", s$run$from), opt("to", s$run$to))
quit(status = as.integer(any(st$status == "ошибка")))

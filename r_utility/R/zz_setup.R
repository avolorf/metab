# Пакеты и загрузка модулей утилиты.

REQUIRED_PACKAGES <- c("data.table", "yaml", "ggplot2", "patchwork", "pheatmap", "caret", "e1071", "kernlab",
                       "randomForest", "pls", "foreach", "doParallel", "jsonlite", "shiny", "rmarkdown", "knitr")
OPTIONAL_PACKAGES <- c(xcms = "Bioconductor: обработка mzML/mzXML/CDF (модуль 1)", MSnbase = "Bioconductor: чтение сырых данных (модуль 1)",
                       IPO = "Bioconductor: оптимизация параметров xcms (модуль 1)", missForest = "CRAN: альтернативная реализация заполнения пропусков")

missing_packages <- function() REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace, TRUE, quietly = TRUE)]

# Первая доступная для записи библиотека пакетов. Если R установлен в Program Files без прав администратора,
# системная библиотека закрыта для записи — тогда создаётся личная (R_LIBS_USER), как предлагает R в интерактивном режиме.
writable_lib <- function() {
  for (l in .libPaths()) if (file.access(l, 2) == 0) return(l)
  l <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(l)) l <- file.path(path.expand("~"), "R", "library")
  l <- strsplit(l, .Platform$path.sep)[[1]][1]
  dir.create(l, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(l, .libPaths()))
  l
}

install_missing <- function(repos = "https://cloud.r-project.org") {
  miss <- missing_packages()
  if (length(miss)) {
    lib <- writable_lib()
    message("Первый запуск: устанавливаются пакеты R (нужен интернет, 5–15 минут): ", paste(miss, collapse = ", "), "\nПапка пакетов: ", lib)
    utils::install.packages(miss, lib = lib, repos = repos)
  }
  miss <- missing_packages()
  if (length(miss)) stop("Не удалось установить пакеты: ", paste(miss, collapse = ", "),
                         ". Проверьте подключение к интернету (или прокси-сервер организации) и запустите программу ещё раз.")
  invisible(TRUE)
}

load_urometabr <- function(root) {
  if (!isTRUE(l10n_info()$`UTF-8`))
    for (loc in c("C.UTF-8", "en_US.UTF-8", "ru_RU.UTF-8", "Russian_Russia.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  options(urometabr.root = normalizePath(root, winslash = "/"), urometabr.dir = file.path(normalizePath(root, winslash = "/"), "R"))
  for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)))
    source(f, encoding = "UTF-8", local = globalenv())   # файлы в UTF-8 (в R 4.0–4.1 под Windows перекодируются)
  invisible(TRUE)
}

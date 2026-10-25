# Пакеты и загрузка модулей утилиты.

REQUIRED_PACKAGES <- c("data.table", "yaml", "ggplot2", "patchwork", "pheatmap", "caret", "e1071", "kernlab",
                       "randomForest", "pls", "foreach", "doParallel", "jsonlite", "shiny", "rmarkdown", "knitr")
OPTIONAL_PACKAGES <- c(xcms = "Bioconductor: обработка mzML/mzXML/CDF (модуль 1)", MSnbase = "Bioconductor: чтение сырых данных (модуль 1)",
                       IPO = "Bioconductor: оптимизация параметров xcms (модуль 1)", missForest = "CRAN: альтернативная реализация заполнения пропусков")

missing_packages <- function() REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace, TRUE, quietly = TRUE)]

install_missing <- function(repos = "https://cloud.r-project.org") {
  miss <- missing_packages()
  if (length(miss)) { message("Установка пакетов: ", paste(miss, collapse = ", ")); utils::install.packages(miss, repos = repos) }
  miss <- missing_packages()
  if (length(miss)) stop("Не удалось установить: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

load_urometabr <- function(root) {
  if (!isTRUE(l10n_info()$`UTF-8`))
    for (loc in c("C.UTF-8", "en_US.UTF-8", "ru_RU.UTF-8", "Russian_Russia.utf8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  options(urometabr.root = normalizePath(root, winslash = "/"), urometabr.dir = file.path(normalizePath(root, winslash = "/"), "R"))
  for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) sys.source(f, envir = globalenv(), keep.source = FALSE)
  invisible(TRUE)
}

# Модуль 0. Загрузка и проверка входных данных: таблица CSV/TSV/TXT (или экспорт MS-DIAL) либо сырые
# mzML/mzXML/CDF, файл метаданных (образец, группа, порядок ввода, повтор). Результат — проверенный
# набор данных или понятное сообщение об ошибке.

read_metadata <- function(path) {
  if (!nzchar(path) || !file.exists(path)) stop("Не найден файл метаданных: ", path)
  md <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8", colClasses = "character")
  names(md) <- tolower(trimws(sub("^﻿", "", names(md))))
  if (!"file" %in% names(md) && "sample" %in% names(md)) md$file <- md$sample
  miss <- setdiff(c("file", "group"), names(md))
  if (length(miss)) stop("В файле метаданных нет обязательных столбцов: ", paste(miss, collapse = ", "),
                         ". Нужны: file (имя файла или строки таблицы), group; необязательно: sample, injection_order, replicate, include, sex, age")
  md$file <- trimws(md$file); md$group <- trimws(md$group)
  if (any(duplicated(md$file))) stop("В метаданных повторяются имена файлов: ", paste(unique(md$file[duplicated(md$file)]), collapse = ", "))
  if (any(md$group == "")) stop("В метаданных не указана группа для: ", paste(md$file[md$group == ""], collapse = ", "))
  md$sample <- if ("sample" %in% names(md)) ifelse(nzchar(md$sample), md$sample, md$file) else md$file
  md$injection_order <- if ("injection_order" %in% names(md)) suppressWarnings(as.numeric(md$injection_order)) else NA_real_
  md$replicate <- if ("replicate" %in% names(md)) suppressWarnings(as.integer(md$replicate)) else 1L
  md$include <- if ("include" %in% names(md)) ifelse(md$include %in% c("0", "FALSE", "false", "нет"), 0L, 1L) else 1L
  if ("age" %in% names(md)) md$age <- suppressWarnings(as.numeric(md$age))
  if ("sex" %in% names(md)) { md$sex <- toupper(trimws(md$sex)); md$sex[md$sex == "М"] <- "M"; md$sex[md$sex == "Ж"] <- "F" }
  md$kind <- ifelse(toupper(md$group) == "QC", "QC", ifelse(grepl("^(blank|бланк)", md$group, ignore.case = TRUE), "blank", "sample"))
  rownames(md) <- md$file
  md
}

# Имя признака → m/z, RT: «RT_mz_аннотация», «mz / RT» (xcms, RT в секундах), «M150.05T2.23»
parse_feature_names <- function(cols) {
  rt <- mz <- rep(NA_real_, length(cols)); nm <- rep("", length(cols))
  m1 <- regmatches(cols, regexec("^([0-9]+\\.?[0-9]*)_([0-9]+\\.?[0-9]*)_?(.*)$", cols))
  m2 <- regmatches(cols, regexec("^([0-9]+\\.?[0-9]*) ?/ ?([0-9]+\\.?[0-9]*)$", cols))
  m3 <- regmatches(cols, regexec("^M([0-9]+\\.?[0-9]*)T([0-9]+\\.?[0-9]*)", cols))
  for (i in seq_along(cols)) {
    if (length(m1[[i]]) == 4) { rt[i] <- as.numeric(m1[[i]][2]); mz[i] <- as.numeric(m1[[i]][3]); nm[i] <- m1[[i]][4] }
    else if (length(m2[[i]]) == 3) { mz[i] <- as.numeric(m2[[i]][2]); rt[i] <- as.numeric(m2[[i]][3]) / 60 }
    else if (length(m3[[i]]) == 3) { mz[i] <- as.numeric(m3[[i]][2]); rt[i] <- as.numeric(m3[[i]][3]) }
  }
  data.frame(id = sprintf("F%05d", seq_along(cols)), original = cols, mz = mz, rt_min = rt,
             name = sub("^low score:\\s*", "", nm), stringsAsFactors = FALSE)
}

read_msdial_table <- function(path) {
  lines <- readLines(path, n = 60, encoding = "UTF-8", warn = FALSE)
  hdr <- which(grepl("^\"?Alignment ID", lines))[1]
  sep <- if (grepl("\t", lines[hdr])) "\t" else ","
  top <- data.table::fread(path, sep = sep, nrows = hdr - 1, header = FALSE, data.table = FALSE, fill = TRUE, encoding = "UTF-8")
  tab <- data.table::fread(path, sep = sep, skip = hdr - 1, header = TRUE, data.table = FALSE, check.names = FALSE, encoding = "UTF-8")
  ft <- trimws(as.character(unlist(top[which(apply(top, 1, function(x) any(trimws(x) == "File type")))[1], ])))
  sc <- which(tolower(ft) %in% c("sample", "qc", "blank", "standard")); sc <- sc[sc <= ncol(tab)]
  X <- t(as.matrix(tab[, sc, drop = FALSE])); storage.mode(X) <- "double"
  info <- data.frame(id = sprintf("F%05d", seq_len(nrow(tab))), original = as.character(tab[["Alignment ID"]]),
                     mz = as.numeric(tab[["Average Mz"]]), rt_min = as.numeric(tab[["Average Rt(min)"]]),
                     name = sub("^low score:\\s*", "", as.character(tab[["Metabolite name"]] %||% "")), stringsAsFactors = FALSE)
  list(X = X, info = info)
}

# Таблица признаков: строки — инъекции (или признаки — определяется по совпадению с метаданными)
read_feature_table <- function(path, md, format = "auto", id_column = "") {
  if (!nzchar(path) || !file.exists(path)) stop("Не найден файл таблицы: ", path)
  head1 <- readLines(path, n = 60, encoding = "UTF-8", warn = FALSE)
  if (format == "msdial" || (format == "auto" && any(grepl("^\"?Alignment ID", head1)))) {
    r <- read_msdial_table(path); r$format <- "MS-DIAL"; return(r)
  }
  raw <- data.table::fread(path, data.table = FALSE, check.names = FALSE, encoding = "UTF-8", colClasses = list(character = 1))
  names(raw)[1] <- sub("^﻿", "", names(raw)[1])
  idc <- if (nzchar(id_column)) id_column else names(raw)[1]
  if (!idc %in% names(raw)) stop("В таблице нет столбца идентификатора '", idc, "'")
  ids <- as.character(raw[[idc]])
  if (sum(ids %in% md$file) == 0 && sum(names(raw) %in% md$file) > 0) {        # признаки по строкам
    sc <- names(raw)[names(raw) %in% md$file]
    X <- t(as.matrix(raw[, sc])); storage.mode(X) <- "double"
    info <- parse_feature_names(ids); return(list(X = X, info = info, format = "признаки по строкам"))
  }
  num <- vapply(raw, function(v) is.numeric(v) || all(is.na(v)), TRUE); num[names(raw) == idc] <- FALSE
  skipped <- setdiff(names(raw)[!num], idc)
  if (length(skipped)) log_msg("Нечисловые столбцы пропущены: ", paste(utils::head(skipped, 10), collapse = ", "))
  X <- as.matrix(raw[, num, drop = FALSE]); storage.mode(X) <- "double"
  rownames(X) <- ids
  dup <- duplicated(ids)
  if (any(dup)) { log_msg("Повторяющиеся строки удалены: ", paste(unique(ids[dup]), collapse = ", ")); X <- X[!dup, , drop = FALSE] }
  list(X = X, info = parse_feature_names(colnames(X)), format = "пробы по строкам")
}

check_raw_file <- function(f) {
  con <- file(f, "rb"); on.exit(close(con)); b <- readBin(con, "raw", 1e6)
  ext <- tolower(tools::file_ext(f)); txt <- rawToChar(b[b != as.raw(0)])
  if (ext == "mzml") {
    if (!grepl("<(indexedmzML|mzML)", txt)) return(list(ok = FALSE, msg = "нет тега <mzML>"))
    n <- regmatches(txt, regexpr("spectrumList count=\"[0-9]+\"", txt))
    list(ok = TRUE, msg = if (length(n)) sub(".*\"([0-9]+)\"", "спектров: \\1", n) else "mzML")
  } else if (ext == "mzxml") {
    if (!grepl("<mzXML", txt)) return(list(ok = FALSE, msg = "нет тега <mzXML>"))
    n <- regmatches(txt, regexpr("scanCount=\"[0-9]+\"", txt))
    list(ok = TRUE, msg = if (length(n)) sub(".*\"([0-9]+)\"", "сканов: \\1", n) else "mzXML")
  } else if (ext %in% c("cdf", "nc")) {
    list(ok = identical(b[1:3], charToRaw("CDF")), msg = "netCDF")
  } else list(ok = FALSE, msg = "неизвестный формат")
}

module0_input <- function(s) {
  out <- module_dir(s, 0)
  md <- read_metadata(project_path(s, s$input$metadata_file))
  chk <- data.frame(check = character(0), status = character(0), detail = character(0))
  add <- function(check, ok, detail) chk[nrow(chk) + 1, ] <<- list(check, if (ok) "ок" else "ошибка", detail)
  add("Метаданные", TRUE, sprintf("%d строк; группы: %s", nrow(md), paste(sprintf("%s %d", names(table(md$group)), table(md$group)), collapse = ", ")))
  if (s$input$type == "raw") {
    d <- project_path(s, s$input$raw_dir)
    if (!dir.exists(d)) stop("Не найдена папка с сырыми данными: ", d)
    files <- list.files(d, pattern = "\\.(mzML|mzXML|CDF)$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
    if (!length(files)) stop("В папке ", d, " нет файлов mzML, mzXML или CDF")
    fc <- lapply(files, check_raw_file)
    bad <- !vapply(fc, `[[`, TRUE, "ok")
    add("Файлы сырых данных", !any(bad), sprintf("%d файлов (%s)%s", length(files), paste(names(table(tolower(tools::file_ext(files)))), collapse = ", "),
                                                 if (any(bad)) paste0("; повреждены: ", paste(basename(files[bad]), collapse = ", ")) else ""))
    if (any(bad)) stop("Повреждённые или неподдерживаемые файлы: ", paste(basename(files[bad]), collapse = ", "))
    key <- tools::file_path_sans_ext(basename(files))
    names(files) <- key
    md$file <- tools::file_path_sans_ext(md$file); rownames(md) <- md$file
    ids <- key
    obj <- list(files = files, meta = NULL, type = "raw")
  } else {
    tb <- read_feature_table(project_path(s, s$input$table_file), md, s$input$table_format, s$input$id_column)
    X <- tb$X; X[!is.finite(X) | X <= 0] <- NA
    colnames(X) <- tb$info$id
    add("Таблица признаков", TRUE, sprintf("%s: %d инъекций × %d признаков; пропусков %.1f%%", tb$format, nrow(X), ncol(X), 100 * mean(is.na(X))))
    if (all(is.na(tb$info$mz))) add("m/z и RT признаков", FALSE, "не распознаны в именах столбцов (ожидается RT_mz_имя, mz / RT или M…T…); анализ продолжится без них")
    ids <- rownames(X)
    obj <- list(X = X, info = tb$info, meta = NULL, type = "table")
  }
  nodata <- setdiff(md$file[md$include == 1], ids); nometa <- setdiff(ids, md$file)
  if (length(nometa)) add("Данные без метаданных", FALSE, paste("исключены:", paste(utils::head(nometa, 20), collapse = ", ")))
  if (length(nodata)) add("Метаданные без данных", FALSE, paste(utils::head(nodata, 20), collapse = ", "))
  use <- intersect(ids, md$file[md$include == 1])
  if (length(use) < 6) stop("Слишком мало совпадений имён между данными и метаданными (", length(use), "). Проверьте столбец file")
  meta <- md[use, , drop = FALSE]
  if (all(is.na(meta$injection_order))) {
    meta$injection_order <- match(use, ids); add("Порядок ввода", FALSE, "не указан — принят порядок строк/файлов")
  } else add("Порядок ввода", !anyNA(meta$injection_order), if (anyNA(meta$injection_order)) "не указан у части файлов" else "указан")
  grp <- table(meta$group[meta$kind == "sample"])
  add("Группы сравнения", length(grp) >= 2 && all(grp >= 3), paste(sprintf("%s %d", names(grp), grp), collapse = ", "))
  if (length(grp) < 2) stop("Нужны минимум две группы образцов (кроме QC и бланков)")
  nqc <- sum(meta$kind == "QC")
  add("Образцы QC", nqc >= 3, sprintf("%d; %s", nqc, if (nqc >= 3) "фильтр RSD будет применён" else "фильтр RSD и детекции в QC будут пропущены"))
  add("Повторные измерения", TRUE, sprintf("%d образцов имеют > 1 инъекции", sum(table(meta$sample[meta$kind == "sample"]) > 1)))
  obj$meta <- meta
  if (obj$type == "table") {
    obj$X <- obj$X[use, , drop = FALSE]
    pct <- 100 * mean(is.na(obj$X[meta$kind == "sample", ]))
    add("Доля пропусков в образцах", pct <= 30, sprintf("%.1f%% (допустимо по ТЗ до 20–30%%)", pct))
  } else obj$files <- obj$files[use]
  write_csv(chk, file.path(out, "input_check.csv"))
  write_csv(meta, file.path(out, "metadata_checked.csv"))
  for (i in seq_len(nrow(chk))) log_msg(sprintf("  [%s] %s: %s", chk$status[i], chk$check[i], chk$detail[i]))
  obj$checks <- chk
  save_module(s, 0, obj)
  if (obj$type == "table") {                                    # для таблиц модуль 1 не нужен: таблица пиков уже есть
    write_csv(cbind(file = rownames(obj$X), as.data.frame(obj$X)), file.path(module_dir(s, 1), "peak_table.csv"))
    write_csv(obj$info, file.path(module_dir(s, 1), "features.csv"))
    save_module(s, 1, obj[c("X", "info", "meta", "type")])
  }
  invisible(obj)
}

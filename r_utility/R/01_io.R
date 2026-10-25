# Загрузка таблицы признаков (CSV/TSV «образцы × признаки», экспорт выравнивания MS-DIAL)
# и метаданных; порядок инъекций.
#
# Результат read_feature_table(): список
#   X     — матрица интенсивностей (строки — инъекции, столбцы — признаки), 0 → NA;
#   info  — data.frame признаков: feature, rt_min, mz, annotation, annotation_conf, ms2;
#   kind  — тип инъекции: sample | QC | blank;
#   inj_order — порядок инъекций, если он есть в файле (MS-DIAL), иначе NULL.

parse_feature_names <- function(cols) {
  # формат «RT_mz_аннотация» (экспорт расчётной таблицы) или «mz / RT» (таблица xcms, RT в секундах)
  m1 <- regmatches(cols, regexec("^([0-9]+\\.?[0-9]*)_([0-9]+\\.?[0-9]*)_(.*)$", cols))
  m2 <- regmatches(cols, regexec("^([0-9]+\\.?[0-9]*) / ([0-9]+\\.?[0-9]*)$", cols))
  rt <- mz <- rep(NA_real_, length(cols)); ann <- rep("Unknown", length(cols))
  for (i in seq_along(cols)) {
    if (length(m1[[i]]) == 4) {
      rt[i] <- as.numeric(m1[[i]][2]); mz[i] <- as.numeric(m1[[i]][3]); ann[i] <- m1[[i]][4]
    } else if (length(m2[[i]]) == 3) {
      mz[i] <- as.numeric(m2[[i]][2]); rt[i] <- as.numeric(m2[[i]][3]) / 60
    }
  }
  feature_info(cols, rt, mz, ann)
}

feature_info <- function(feature, rt, mz, ann, ms2 = NA) {
  low <- grepl("^low score:", ann)
  clean <- sub("^low score:\\s*", "", ann)
  conf <- ifelse(clean %in% c("Unknown", "unknown", ""), "unknown", ifelse(low, "low", "high"))
  data.frame(feature = feature, rt_min = rt, mz = mz, annotation = clean, annotation_conf = conf,
             ms2 = ms2, stringsAsFactors = FALSE, row.names = feature)
}

classify_injections <- function(ids, groups, cfg) {
  kind <- rep("sample", length(ids))
  kind[grepl(cfg$design$qc_pattern, ids, ignore.case = TRUE)] <- "QC"
  bl <- grepl(cfg$design$blank_pattern, ids, ignore.case = TRUE, perl = TRUE) |
        grepl(cfg$design$blank_pattern, groups, ignore.case = TRUE, perl = TRUE)
  kind[bl] <- "blank"
  stats::setNames(kind, ids)
}

read_wide <- function(path, cfg) {
  raw <- data.table::fread(path, encoding = "UTF-8", colClasses = list(character = cfg$input$id_column),
                           check.names = FALSE, data.table = FALSE)
  id <- cfg$input$id_column
  if (!id %in% names(raw)) stop("В таблице нет столбца идентификатора проб '", id, "'")
  mcols <- intersect(cfg$input$meta_columns, names(raw))
  fcols <- setdiff(names(raw), mcols)
  ndup <- sum(duplicated(raw[[id]]))
  raw <- raw[!duplicated(raw[[id]]), , drop = FALSE]
  X <- as.matrix(raw[, fcols, drop = FALSE]); storage.mode(X) <- "double"
  X[X == 0] <- NA
  rownames(X) <- raw[[id]]
  grp <- if ("Group" %in% names(raw)) as.character(raw$Group) else rep("", nrow(raw))
  list(X = X, info = parse_feature_names(fcols), kind = classify_injections(rownames(X), grp, cfg),
       inj_order = NULL, n_dup = ndup)
}

read_msdial <- function(path, cfg) {
  # Экспорт MS-DIAL «Alignment result → Area/Height»: над заголовком строки Class, File type,
  # Injection order, Batch ID; признаки — по строкам.
  lines <- readLines(path, n = 50, encoding = "UTF-8")
  sep <- if (grepl("\t", lines[1])) "\t" else ","
  hdr <- which(grepl("^\"?Alignment ID", lines))[1]
  if (is.na(hdr)) stop("Не найден заголовок 'Alignment ID' — это не экспорт MS-DIAL")
  top <- data.table::fread(path, sep = sep, nrows = hdr - 1, header = FALSE, data.table = FALSE,
                           encoding = "UTF-8", fill = TRUE)
  tab <- data.table::fread(path, sep = sep, skip = hdr - 1, header = TRUE, data.table = FALSE,
                           encoding = "UTF-8", check.names = FALSE)
  lab <- function(key) { r <- which(apply(top, 1, function(x) any(trimws(x) == key))); if (length(r)) top[r[1], ] else NULL }
  ftype <- lab("File type"); inj <- lab("Injection order"); cls <- lab("Class")
  if (is.null(ftype)) stop("В экспорте MS-DIAL нет строки 'File type'")
  ft <- trimws(as.character(unlist(ftype)))
  scols <- which(tolower(ft) %in% c("sample", "qc", "blank", "standard"))
  scols <- scols[scols <= ncol(tab)]
  rt <- as.numeric(tab[["Average Rt(min)"]]); mz <- as.numeric(tab[["Average Mz"]])
  ann <- as.character(tab[["Metabolite name"]] %||% rep("Unknown", nrow(tab)))
  ms2 <- if ("MS/MS assigned" %in% names(tab)) as.logical(toupper(tab[["MS/MS assigned"]])) else NA
  fid <- make.unique(sprintf("%.3f_%.5f_%s", rt, mz, ann))
  X <- t(as.matrix(tab[, scols, drop = FALSE])); storage.mode(X) <- "double"
  X[X == 0] <- NA
  colnames(X) <- fid
  ids <- names(tab)[scols]; rownames(X) <- ids
  k <- tolower(ft[scols])
  kind <- stats::setNames(ifelse(k == "qc", "QC", ifelse(k == "blank", "blank", "sample")), ids)
  pat <- classify_injections(ids, rep("", length(ids)), cfg)     # имя пробы важнее пометки File type
  over <- kind == "sample" & pat != "sample"
  kind[over] <- pat[over]
  io <- if (!is.null(inj)) stats::setNames(suppressWarnings(as.numeric(unlist(inj)[scols])), ids) else NULL
  list(X = X, info = feature_info(fid, rt, mz, ann, ms2), kind = kind, inj_order = io, n_dup = 0L,
       msdial_class = if (!is.null(cls)) stats::setNames(as.character(unlist(cls)[scols]), ids) else NULL)
}

read_feature_table <- function(cfg) {
  f <- cfg$input$format
  if (f == "wide") return(read_wide(cfg$input$data, cfg))
  if (f == "msdial") return(read_msdial(cfg$input$data, cfg))
  if (f == "xcms") return(xcms_peak_table(cfg))      # R/01b_xcms.R (нужен Bioconductor-пакет xcms)
  stop("Неизвестный формат: ", f)
}

load_metadata <- function(path, kind, cfg) {
  md <- data.table::fread(path, colClasses = list(character = "sample"), data.table = FALSE, encoding = "UTF-8")
  need <- c("sample", "group", "sex", "age", "include")
  miss <- setdiff(need, names(md))
  if (length(miss)) stop("В ", path, " нет столбцов: ", paste(miss, collapse = ", "))
  md$sex <- toupper(trimws(ifelse(is.na(md$sex), "", md$sex)))
  md$sex[md$sex == "М"] <- "M"; md$sex[md$sex == "Ж"] <- "F"
  md$age <- suppressWarnings(as.numeric(md$age))
  md$include <- ifelse(is.na(md$include), 0L, as.integer(md$include))
  rownames(md) <- md$sample
  ids <- names(kind)
  smp <- ids[kind == "sample"]
  if (length(setdiff(smp, md$sample)))
    log_msg("ВНИМАНИЕ: пробы без метаданных (исключены): ", paste(setdiff(smp, md$sample), collapse = ", "))
  nd <- setdiff(md$sample[md$include == 1], smp)
  if (length(nd)) log_msg("ВНИМАНИЕ: в метаданных, но нет в данных: ", paste(nd, collapse = ", "))
  meta <- data.frame(sample = ids, kind = unname(kind), stringsAsFactors = FALSE, row.names = ids)
  for (c in c("group", "sex", "age", "include", "exclude_reason", "injection_order"))
    if (c %in% names(md)) meta[[c]] <- md[ids, c]
  meta$group[meta$kind == "QC"] <- "QC"
  meta$group[meta$kind == "blank"] <- "Blank"
  meta$include[is.na(meta$include) | meta$kind != "sample"] <- 0L
  bad <- meta$include == 1 & !meta$group %in% cfg$design$groups
  if (any(bad)) {
    log_msg("ВНИМАНИЕ: include=1, но группа вне design$groups (исключены): ", paste(ids[bad], collapse = ", "))
    meta$include[bad] <- 0L
  }
  meta$study <- meta$kind == "sample" & meta$include == 1
  meta
}

injection_order <- function(meta, file_order = NULL, qc_every = 20) {
  # 1) столбец injection_order в метаданных; 2) порядок из MS-DIAL; 3) правило: QC1 первым,
  # затем QC после каждых qc_every образцов (образцы по номеру).
  if (!is.null(meta$injection_order) && any(!is.na(meta$injection_order)) && all(!is.na(meta$injection_order[meta$kind != "blank"])))
    return(stats::setNames(as.numeric(meta$injection_order), meta$sample))
  if (!is.null(file_order) && all(!is.na(file_order))) return(file_order[meta$sample])
  num <- function(s) suppressWarnings(as.numeric(gsub("\\D", "", s)))
  qcs <- meta$sample[meta$kind == "QC"]; qcs <- qcs[order(num(qcs))]
  smp <- meta$sample[meta$kind == "sample"]; smp <- smp[order(nchar(smp), smp)]
  ord <- stats::setNames(rep(NA_real_, nrow(meta)), meta$sample)
  pos <- 0; qi <- 1
  if (length(qcs)) { ord[qcs[1]] <- pos; qi <- 2 }
  for (i in seq_along(smp)) {
    pos <- pos + 1; ord[smp[i]] <- pos
    if (i %% qc_every == 0 && qi <= length(qcs)) { pos <- pos + 1; ord[qcs[qi]] <- pos; qi <- qi + 1 }
  }
  ord
}

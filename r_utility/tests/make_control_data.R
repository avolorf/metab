# Генератор контрольных наборов данных с известной «истиной» (ground truth).
# Моделируются: группы Control/BC/PC с дисбалансом по возрасту и полу, разведение мочи, технический шум,
# пропуски ниже предела обнаружения, нестабильные и фоновые (бланковые) признаки, аддукты одного
# соединения, признаки, зависящие только от возраста, ряд экзогенных диэтаноламинов (чаще в PC),
# «лекарственный» признак, повышенный в BC, и истинные отличия BC от контроля.
#
#   Rscript tests/make_control_data.R <каталог> [seed] [null]
# null = TRUE: все групповые эффекты нулевые (проверка ложноположительных результатов).

make_control_data <- function(dir, seed = 1, null = FALSE, series_file = NULL, n_feat = 1500,
                              n = c(Control = 45, BC = 29, PC = 20), n_qc = 8, n_blank = 4) {
  set.seed(seed)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ELEM <- c(C = 12, H = 1.00782503, N = 14.0030740, O = 15.9949146)
  fm <- function(f) { m <- regmatches(f, gregexpr("([A-Z][a-z]?)(\\d*)", f))[[1]]
    sum(vapply(m, function(t) ELEM[[sub("\\d+$", "", t)]] * max(1, as.numeric(sub("^[A-Za-z]+", "", t)), na.rm = TRUE), 0)) }
  # --- пробы
  grp <- rep(names(n), n); ns <- length(grp)
  sex <- ifelse(grp == "PC", "M", ifelse(grp == "BC", ifelse(runif(ns) < 0.83, "M", "F"), ifelse(runif(ns) < 0.42, "M", "F")))
  age <- round(ifelse(grp == "Control", rnorm(ns, 53, 14), ifelse(grp == "BC", rnorm(ns, 63, 10), rnorm(ns, 70, 7))))
  age <- pmin(pmax(age, 18), 88)
  ids <- as.character(seq_len(ns))[sample(ns)]         # номера проб не связаны с группой
  dil <- exp(rnorm(ns, 0, 0.5))                         # истинный фактор разведения
  contam <- pmax(0, rnorm(ns, ifelse(grp == "PC", 1.5, ifelse(grp == "BC", 0.6, 0.2)), 0.7))
  # --- признаки
  nf <- n_feat
  type <- rep("stable", nf)
  type[1:30] <- "BC_up"; type[31:50] <- "BC_down"; type[51:70] <- "age_only"; type[71:170] <- "noisy"
  type[171:270] <- "rare"; type[271:320] <- "blank"; type[321:323] <- "drug_BC"
  rt <- round(runif(nf, 0.8, 20), 3); mz <- round(runif(nf, 150, 900), 5)
  base <- rnorm(nf, 15, 1.8)
  tech <- ifelse(type == "noisy", 0.75, 0.08)           # SD log2 технического шума (noisy → RSD > 30 %)
  eff_bc <- ifelse(type == "BC_up", 0.7, ifelse(type == "BC_down", -0.7, ifelse(type == "drug_BC", 2, 0)))
  if (null) eff_bc[] <- 0
  bio <- matrix(rnorm(ns * nf, 0, 0.45), ns, nf)
  Ltrue <- matrix(base, ns, nf, byrow = TRUE) + bio + outer(grp == "BC", eff_bc) +
           outer((age - 55) * 0.04, as.numeric(type == "age_only"))
  # сырые значения = истинные × разведение × технический шум
  raw <- 2^(Ltrue + matrix(rnorm(ns * nf, 0, rep(tech, each = ns)), ns, nf)) * dil
  # аддукты: признаки 400–499 получают Na-аддукт (+21.98194, тот же RT, интенсивность ×0.3)
  src <- 400:499
  add <- raw[, src] * 0.3 * 2^matrix(rnorm(ns * length(src), 0, 0.05), ns)
  # ряд экзогенных соединений
  ser <- utils::read.csv(series_file, stringsAsFactors = FALSE)
  ser <- ser[ser$in_index == 1, ]
  sraw <- sapply(seq_len(nrow(ser)), function(k) 2^(14 + 1.2 * contam + rnorm(ns, 0, 0.3)) * dil)
  X <- cbind(raw, add, sraw)
  fmz <- c(mz, mz[src] + 21.981942, vapply(ser$formula, fm, 0) + 1.007276)
  frt <- c(rt, rt[src] + 0.001, ser$rt_min)
  ann <- c(rep("Unknown", nf), rep("Unknown", length(src)), rep("Unknown", nrow(ser)))
  ann[321:323] <- c("Acetaminophen", "low score: Paracetamol sulfate", "4-Acetamidophenol glucuronide")
  ann[sample(500:nf, 60)] <- paste("low score: putative compound", 1:60)
  ftype <- c(type, rep("adduct", length(src)), rep("series", nrow(ser)))
  # пропуски: ниже предела обнаружения и «редкие» признаки
  lod <- 2^13
  X[X < lod] <- NA
  rare <- which(ftype == "rare"); X[, rare][runif(ns * length(rare)) < 0.75] <- NA
  # QC = средний профиль проб + малый шум; бланки: фон только у признаков типа blank
  qc <- t(replicate(n_qc, colMeans(X, na.rm = TRUE) * 2^rnorm(ncol(X), 0, 0.06)))
  qc[, rare] <- ifelse(runif(n_qc * length(rare)) < 0.5, NA, qc[, rare])
  qc[, ftype == "noisy"] <- qc[, ftype == "noisy"] * 2^matrix(rnorm(n_qc * sum(ftype == "noisy"), 0, 0.75), n_qc)
  bl <- matrix(NA_real_, n_blank, ncol(X))
  bl[, ftype == "blank"] <- qc[1:n_blank, ftype == "blank"] * 0.6
  # таблица в формате «RT_mz_аннотация»
  fnames <- make.unique(sprintf("%.3f_%.5f_%s", frt, fmz, ann))
  qids <- paste0("QC", seq_len(n_qc)); bids <- paste0("Blank", seq_len(n_blank))
  M <- rbind(X, qc, bl); M[is.na(M)] <- 0
  tab <- data.frame(Sample = c(ids, qids, bids), Group = c(grp, rep("QC", n_qc), rep("Blank", n_blank)),
                    signif(M, 6), check.names = FALSE, stringsAsFactors = FALSE)
  names(tab)[-(1:2)] <- fnames
  # порядок ввода: QC в начале и через равные интервалы, бланки в начале и в конце
  step <- ceiling(ns / (n_qc - 1))
  seqv <- c(bids[1], qids[1]); qi <- 2
  perm <- sample(ns)                                    # рандомизированный порядок ввода
  for (i in seq_len(ns)) { seqv <- c(seqv, ids[perm[i]]); if (i %% step == 0 && qi <= n_qc) { seqv <- c(seqv, qids[qi]); qi <- qi + 1 } }
  seqv <- c(seqv, qids[qi:n_qc][qi <= n_qc], bids[-1])
  seqv <- unique(seqv)
  md <- data.frame(sample = c(ids, qids, bids), group = c(grp, rep("QC", n_qc), rep("Blank", n_blank)),
                   sex = c(sex, rep("", n_qc + n_blank)), age = c(age, rep(NA, n_qc + n_blank)),
                   include = c(rep(1, ns), rep(0, n_qc + n_blank)),
                   injection_order = match(c(ids, qids, bids), seqv), stringsAsFactors = FALSE)
  data.table::fwrite(tab, file.path(dir, "feature_table.csv"))
  data.table::fwrite(md, file.path(dir, "metadata.csv"))
  truth <- data.frame(feature = fnames, type = ftype, effect_BC = c(eff_bc, rep(0, length(src) + nrow(ser))),
                      parent = c(rep(NA, nf), fnames[src], rep(NA, nrow(ser))), stringsAsFactors = FALSE)
  data.table::fwrite(truth, file.path(dir, "truth.csv"))
  data.table::fwrite(data.frame(sample = ids, dilution_true = dil, contam_true = contam), file.path(dir, "truth_samples.csv"))
  invisible(list(table = tab, meta = md, truth = truth))
}

# Тот же набор в формате экспорта MS-DIAL (признаки по строкам, 4 строки заголовка над таблицей)
write_msdial_export <- function(dir) {
  tab <- data.table::fread(file.path(dir, "feature_table.csv"), data.table = FALSE, check.names = FALSE, colClasses = list(character = "Sample"))
  md <- data.table::fread(file.path(dir, "metadata.csv"), data.table = FALSE, colClasses = list(character = "sample"))
  f <- names(tab)[-(1:2)]
  parts <- do.call(rbind, regmatches(f, regexec("^([0-9.]+)_([0-9.]+)_(.*)$", f)))
  nm <- sub("\\.[0-9]+$", "", parts[, 4])
  info <- data.frame(`Alignment ID` = seq_along(f) - 1, `Average Rt(min)` = parts[, 2], `Average Mz` = parts[, 3],
                     `Metabolite name` = nm, `Adduct type` = "[M+H]+", `MS/MS assigned` = "False", check.names = FALSE)
  vals <- t(as.matrix(tab[, -(1:2)])); colnames(vals) <- tab$Sample
  ft <- ifelse(grepl("^QC", tab$Sample), "QC", ifelse(grepl("^Blank", tab$Sample), "Blank", "Sample"))
  io <- md$injection_order[match(tab$Sample, md$sample)]
  pad <- rep("", ncol(info) - 1)
  hdr <- rbind(c(pad, "Class", tab$Group), c(pad, "File type", ft), c(pad, "Injection order", io), c(pad, "Batch ID", rep(1, nrow(tab))))
  out <- file.path(dir, "msdial_area.txt")
  con <- file(out, "w", encoding = "UTF-8")
  for (i in seq_len(nrow(hdr))) writeLines(paste(hdr[i, ], collapse = "\t"), con)
  writeLines(paste(c(names(info), tab$Sample), collapse = "\t"), con)
  body <- cbind(as.matrix(info), vals)
  for (i in seq_len(nrow(body))) writeLines(paste(body[i, ], collapse = "\t"), con)
  close(con)
  out
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
  make_control_data(a[1], seed = if (length(a) > 1) as.numeric(a[2]) else 1, null = length(a) > 2 && as.logical(a[3]),
                    series_file = file.path(here, "..", "example", "contaminant_series.csv"))
}

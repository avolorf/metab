# Контрольные наборы данных для испытаний.
#  Набор А* — синтетический заменитель набора этапа 1 (106 образцов: 50 CG, 56 BC; QC; повторные вводы),
#             используется, пока реальный набор А не подставлен (tests/ispytaniya.R --setA=...).
#  Набор В  — из любого набора (А или А*): метки групп перемешиваются (реальные различия разрушаются),
#             в n_spike признаков вносятся различия известной величины, доля пропусков доводится до заданной.
#  Форматы  — та же таблица в CSV, TSV, TXT, «признаки по строкам»; минимальные mzML/mzXML для модуля 0.

new_project <- function(dir, table_file = "", metadata_file = "", type = "table", raw_dir = "", ...) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  s <- default_settings()
  s$project$name <- basename(dir)
  s$input$type <- type; s$input$table_file <- table_file; s$input$metadata_file <- metadata_file; s$input$raw_dir <- raw_dir
  over <- list(...); if (length(over)) s <- merge_settings(s, over)
  f <- file.path(dir, "settings.yaml"); write_settings(s, f); f
}

make_set_A_surrogate <- function(dir, seed = 1, n = c(CG = 50, BC = 56), n_feat = 1200, n_qc = 8, n_rep = 10, n_true = 60) {
  set.seed(seed); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  grp <- rep(names(n), n); ns <- length(grp); smp <- sprintf("S%03d", seq_len(ns))
  base <- rnorm(n_feat, 14, 2); bio <- matrix(rnorm(ns * n_feat, 0, 0.6), ns)
  eff <- numeric(n_feat); tr <- sample(n_feat, n_true); eff[tr] <- sample(c(-1, 1), n_true, TRUE) * runif(n_true, 0.3, 0.8)
  tech <- rep(0.12, n_feat); noisy <- sample(setdiff(seq_len(n_feat), tr), round(0.15 * n_feat)); tech[noisy] <- 0.6
  L <- matrix(base, ns, n_feat, byrow = TRUE) + bio + outer(grp == "BC", eff)
  drift <- outer(seq(-0.15, 0.15, length.out = ns), rep(1, n_feat))               # дрейф чувствительности по порядку ввода
  ord <- sample(ns)
  inj <- function(Lm, k) 2^(Lm + matrix(rnorm(length(Lm), 0, rep(tech, each = nrow(Lm))), nrow(Lm)) + drift[k, , drop = FALSE])
  X <- inj(L[ord, ], seq_len(ns)); files <- smp[ord]; groups <- grp[ord]
  rep_i <- sample(seq_len(ns), n_rep)                                               # повторные вводы части образцов
  Xr <- inj(L[ord[rep_i], , drop = FALSE], rep_i)
  qcL <- matrix(colMeans(L), n_qc, n_feat, byrow = TRUE)
  Xq <- 2^(qcL + matrix(rnorm(n_qc * n_feat, 0, rep(tech, each = n_qc)), n_qc))
  M <- rbind(X, Xr, Xq)
  lod <- stats::quantile(M, 0.12); M[M < lod] <- NA                                # пропуски ниже предела обнаружения
  rare <- sample(n_feat, round(0.1 * n_feat)); M[, rare][runif(nrow(M) * length(rare)) < 0.6] <- NA
  mz <- round(runif(n_feat, 150, 1000), 4); rt <- round(runif(n_feat, 30, 1300), 1)
  ids <- c(files, paste0(files[rep_i], "_r2"), sprintf("QC%02d", seq_len(n_qc)))
  tab <- data.frame(file = ids, M, check.names = FALSE); names(tab)[-1] <- sprintf("%.4f / %.1f", mz, rt)
  # порядок ввода: QC в начале и через равные интервалы, повторы — в конце серии
  seqv <- c("QC01", unlist(lapply(split(c(files, paste0(files[rep_i], "_r2")), ceiling(seq_len(ns + n_rep) / ceiling((ns + n_rep) / (n_qc - 1)))),
                                  function(b) b)))
  qpos <- round(seq(1, length(seqv), length.out = n_qc))[-1]
  for (k in rev(seq_along(qpos))) seqv <- append(seqv, sprintf("QC%02d", k + 1), qpos[k])
  md <- data.frame(file = ids, sample = c(files, files[rep_i], ids[(ns + n_rep + 1):length(ids)]),
                   group = c(groups, groups[rep_i], rep("QC", n_qc)), injection_order = match(ids, seqv),
                   replicate = c(rep(1, ns), rep(2, n_rep), rep(1, n_qc)))
  data.table::fwrite(tab, file.path(dir, "peak_table.csv")); data.table::fwrite(md, file.path(dir, "metadata.csv"))
  data.table::fwrite(data.frame(column = names(tab)[-1], true_effect = eff, noisy = seq_len(n_feat) %in% noisy), file.path(dir, "truth_A.csv"))
  invisible(list(table = file.path(dir, "peak_table.csv"), metadata = file.path(dir, "metadata.csv")))
}

make_set_V <- function(table_file, metadata_file, dir, n_spike = 15, effect = c(0.8, 1.5), missing = 0.25, seed = 2,
                       case = NULL, control = NULL) {
  set.seed(seed); dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tab <- data.table::fread(table_file, data.table = FALSE, check.names = FALSE, colClasses = list(character = 1))
  md <- data.table::fread(metadata_file, data.table = FALSE, colClasses = "character")
  names(md) <- tolower(names(md)); if (is.null(md$sample)) md$sample <- md$file
  X <- as.matrix(tab[, -1]); X[X == 0] <- NA; rownames(X) <- tab[[1]]
  isq <- toupper(md$group) == "QC" | grepl("^blank", md$group, ignore.case = TRUE)
  st <- unique(md[!isq, c("sample", "group")])
  st$group <- sample(st$group)                                         # перемешивание меток на уровне образца
  md$group[!isq] <- st$group[match(md$sample[!isq], st$sample)]
  gl <- sort(unique(st$group)); case <- case %||% gl[1]
  rows <- tab[[1]][!isq[match(tab[[1]], md$file)]]
  inS <- rownames(X) %in% md$file[!isq]
  det <- colMeans(!is.na(X[inS, ])); qcr <- rownames(X) %in% md$file[isq]
  rsd <- apply(X[qcr, , drop = FALSE], 2, function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE))
  cand <- which(det >= 0.95 & !is.na(rsd) & rsd <= 0.2)
  sp <- sample(cand, n_spike); e <- round(runif(n_spike, effect[1], effect[2]), 2) * sample(c(-1, 1), n_spike, TRUE)
  rc <- rownames(X) %in% md$file[md$group == case]
  X[rc, sp] <- sweep(X[rc, sp, drop = FALSE], 2, 2^e, "*")
  truth <- data.frame(column = colnames(X)[sp], log2_effect = e, case_group = case)
  # пропуски: доводим долю пропусков в образцах до заданной (половина — у наименьших значений, половина — случайно)
  S <- X[inS, ]; need <- round(missing * length(S)) - sum(is.na(S))
  mask <- matrix(FALSE, nrow(S), ncol(S))
  if (need > 0) {
    obs <- which(!is.na(S)); low <- obs[order(S[obs])][seq_len(round(need / 2))]
    rnd <- sample(setdiff(obs, low), need - length(low)); mask[c(low, rnd)] <- TRUE
  }
  true_vals <- S[mask]; S[mask] <- NA; X[inS, ] <- S
  tab2 <- data.frame(file = rownames(X), X, check.names = FALSE)
  data.table::fwrite(tab2, file.path(dir, "peak_table.csv")); data.table::fwrite(md, file.path(dir, "metadata.csv"))
  data.table::fwrite(truth, file.path(dir, "truth_V.csv"))
  saveRDS(list(mask_rows = rownames(S), mask = mask, true_values = true_vals), file.path(dir, "introduced_missing.rds"))
  invisible(list(table = file.path(dir, "peak_table.csv"), metadata = file.path(dir, "metadata.csv"), truth = truth,
                 missing_pct = 100 * mean(is.na(S))))
}

write_format_variants <- function(table_file, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tab <- data.table::fread(table_file, data.table = FALSE, check.names = FALSE, colClasses = list(character = 1))
  data.table::fwrite(tab, file.path(dir, "table.csv"))
  data.table::fwrite(tab, file.path(dir, "table.tsv"), sep = "\t")
  data.table::fwrite(tab, file.path(dir, "table.txt"), sep = "\t")
  data.table::fwrite(tab, file.path(dir, "table_semicolon.txt"), sep = ";")
  tr <- data.frame(feature = names(tab)[-1], t(as.matrix(tab[, -1])), check.names = FALSE); names(tr)[-1] <- tab[[1]]
  data.table::fwrite(tr, file.path(dir, "table_features_in_rows.csv"))
  invisible(dir)
}

# Минимальные корректные mzML и mzXML (заголовок, один спектр) — для проверки модуля 0
write_raw_stubs <- function(dir, ids) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_along(ids)) {
    if (i %% 2) writeLines(c('<?xml version="1.0" encoding="utf-8"?>', '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">',
      '<mzML version="1.1.0"><run id="r"><spectrumList count="1" defaultDataProcessingRef="dp">',
      '<spectrum index="0" id="scan=1" defaultArrayLength="0"/></spectrumList></run></mzML></indexedmzML>'), file.path(dir, paste0(ids[i], ".mzML")))
    else writeLines(c('<?xml version="1.0" encoding="ISO-8859-1"?>', '<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_3.2">',
      '<msRun scanCount="1"><scan num="1" msLevel="1" peaksCount="0"/></msRun></mzXML>'), file.path(dir, paste0(ids[i], ".mzXML")))
  }
  invisible(dir)
}

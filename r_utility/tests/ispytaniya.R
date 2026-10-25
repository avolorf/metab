#!/usr/bin/env Rscript
# Функциональные испытания утилиты (10 пунктов программы испытаний). Протокол, журналы и sessionInfo()
# каждого прогона собираются в папку --out.
#
#   Rscript tests/ispytaniya.R --out=протокол
#          [--setA=таблица.csv,метаданные.csv]   набор А (этап 1, 106 образцов); без него — синтетический заменитель А*
#          [--setB=таблица.csv,метаданные.csv]   набор Б (когорта 2025–2026, 94 образца)
#          [--raw=папка_mzML,метаданные.csv]     сырые данные для п. 1 (нужен xcms)
#          [--cores=4]
if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows")
  for (loc in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d = NULL) { a <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(a)) sub("^[^=]*=", "", a) else d }
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(normalizePath(self)), ".."))
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8"); load_urometabr(root)
source(file.path(root, "tests", "make_test_data.R"), encoding = "UTF-8")
OUT <- normalizePath(arg("out", "ispytaniya_out"), mustWork = FALSE); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
W <- file.path(OUT, "work"); dir.create(W, showWarnings = FALSE)
CORES <- as.integer(arg("cores", "1"))
P <- data.frame(trial = integer(0), item = character(0), value = character(0), target = character(0), result = character(0))
add <- function(trial, item, value, target = "", result = "") {
  P[nrow(P) + 1, ] <<- list(trial, item, as.character(value), target, result)
  data.table::fwrite(P, file.path(OUT, "protocol.csv"), bom = TRUE)
  cat(sprintf("[%d] %s: %s %s %s\n", trial, item, value, if (nzchar(target)) paste0("(цель: ", target, ")") else "", result))
}
ok <- function(x) if (isTRUE(x)) "пройден" else "не пройден"

# Прогон проекта; журналы и sessionInfo копируются в протокол
run_project <- function(name, table, meta, from = 0, to = 9, ...) {
  d <- file.path(W, name); unlink(d, recursive = TRUE)
  f <- new_project(d, normalizePath(table), normalizePath(meta), project = list(cores = CORES), ...)
  s <- read_settings(f)
  t0 <- proc.time()[["elapsed"]]
  r <- NULL
  tryCatch(utils::capture.output(r <- run_pipeline(s, from, to)), error = function(e) NULL)
  st <- if (!is.null(r)) r else data.frame(module = NA, status = "ошибка", time_s = NA, peak_ram_mb = NA)
  dl <- file.path(OUT, "runs", name); dir.create(dl, showWarnings = FALSE, recursive = TRUE)
  file.copy(list.files(file.path(d, "logs"), full.names = TRUE), dl, overwrite = TRUE)
  list(dir = d, s = s, status = st, ok = all(st$status %in% c("выполнен", "пропущен")), time = proc.time()[["elapsed"]] - t0,
       ram = suppressWarnings(max(st$peak_ram_mb, na.rm = TRUE)))
}
rd <- function(r, id, f) as.data.frame(data.table::fread(file.path(r$dir, MODULES$dir[MODULES$id == id], f), encoding = "UTF-8"))
md5 <- function(r, files) vapply(files, function(f) digest::digest(file = file.path(r$dir, f), algo = "md5"), "")

# ------------------------------------------------------------------ наборы данных
A <- arg("setA"); B <- arg("setB")
if (is.null(A)) { a <- make_set_A_surrogate(file.path(W, "setA_surrogate")); A <- c(a$table, a$metadata); A_real <- FALSE
} else { A <- strsplit(A, ",")[[1]]; A_real <- TRUE }
if (!is.null(B)) B <- strsplit(B, ",")[[1]]
add(0, "Набор А", if (A_real) paste("реальный:", basename(A[1])) else "синтетический заменитель А* (реальный набор не передан)")
add(0, "Набор Б", if (is.null(B)) "не передан — пункты для набора Б не выполнялись" else basename(B[1]))

# ------------------------------------------------------------------ 1. форматы входа
fv <- write_format_variants(A[1], file.path(W, "formats"))
ref <- NULL
for (f in c("table.csv", "table.tsv", "table.txt", "table_semicolon.txt", "table_features_in_rows.csv")) {
  r <- run_project(paste0("fmt_", tools::file_path_sans_ext(f)), file.path(fv, f), A[2], 0, 3)
  X <- if (r$ok) load_module(r$s, 3)$X else NULL
  same <- if (is.null(ref)) { ref <- X; TRUE } else isTRUE(all.equal(unname(X), unname(ref)))
  add(1, paste("Формат", f), if (r$ok) sprintf("модули 0–3 выполнены; %d × %d", nrow(X), ncol(X)) else "ошибка", "без ошибок, одинаковый результат",
      ok(r$ok && same))
}
ids <- data.table::fread(A[2], data.table = FALSE, colClasses = "character")$file
rawd <- file.path(W, "raw_stub"); write_raw_stubs(rawd, ids)
d <- file.path(W, "fmt_raw_stub"); unlink(d, recursive = TRUE)
s <- read_settings(new_project(d, "", normalizePath(A[2]), type = "raw", raw_dir = rawd))
invisible(log_start(s)); m0 <- tryCatch({ invisible(utils::capture.output(module0_input(s))); TRUE }, error = function(e) conditionMessage(e))
add(1, "Формат mzML и mzXML: модуль 0 (проверка файлов)", if (isTRUE(m0)) sprintf("%d файлов распознаны", length(ids)) else m0, "файлы распознаны", ok(isTRUE(m0)))
RAW <- arg("raw")
if (!requireNamespace("xcms", quietly = TRUE)) {
  add(1, "Формат mzML и mzXML: модуль 1 (xcms)", "не выполнен: пакет xcms не установлен в среде испытаний", "таблица пиков", "не выполнен")
} else if (is.null(RAW)) {
  add(1, "Формат mzML и mzXML: модуль 1 (xcms)", "не выполнен: не передан --raw", "таблица пиков", "не выполнен")
} else {
  rr <- strsplit(RAW, ",")[[1]]; d <- file.path(W, "fmt_raw"); unlink(d, recursive = TRUE)
  r <- run_project("fmt_raw", "", rr[2], 0, 3, input = list(type = "raw", raw_dir = normalizePath(rr[1])))
  add(1, "Формат mzML/mzXML: модули 0–3 по сырым данным", if (r$ok) sprintf("выполнено за %.0f с", r$time) else "ошибка", "без ошибок", ok(r$ok))
}

# ------------------------------------------------------------------ 2. объём сессии; 8. время
subset_set <- function(tab, meta, n, name) {
  t <- data.table::fread(tab, data.table = FALSE, check.names = FALSE, colClasses = list(character = 1))
  m <- data.table::fread(meta, data.table = FALSE, colClasses = "character"); names(m) <- tolower(names(m)); if (is.null(m$sample)) m$sample <- m$file
  isq <- toupper(m$group) == "QC" | grepl("^blank", m$group, ignore.case = TRUE)
  st <- unique(m[!isq, c("sample", "group")]); set.seed(n)
  keep <- unlist(lapply(split(st$sample, st$group), function(x) x[seq_len(round(n * length(x) / nrow(st)))]))
  m2 <- m[isq | m$sample %in% keep, ]; t2 <- t[t[[1]] %in% m2$file, ]
  d <- file.path(W, name); dir.create(d, showWarnings = FALSE)
  data.table::fwrite(t2, file.path(d, "peak_table.csv")); data.table::fwrite(m2, file.path(d, "metadata.csv"))
  c(file.path(d, "peak_table.csv"), file.path(d, "metadata.csv"))
}
for (n in c(50, 100)) {
  ss <- subset_set(A[1], A[2], n, paste0("subset", n))
  r <- run_project(paste0("volume_", n), ss[1], ss[2])
  add(2, sprintf("Объём %d образцов (набор А%s)", n, if (A_real) "" else "*"),
      sprintf("%s; время %.1f мин; пиковая память R %.0f МБ", if (r$ok) "выполнено" else "ошибка", r$time / 60, r$ram), "50–100 образцов без ошибок", ok(r$ok))
  data.table::fwrite(r$status, file.path(OUT, "runs", paste0("volume_", n), "modules.csv"), bom = TRUE)
}

# ------------------------------------------------------------------ 3. устойчивые признаки; 7. сходимость; 8. время
rA <- run_project("setA_full", A[1], A[2])
nfeat <- function(r) ncol(load_module(r$s, 3)$X)
nA <- nfeat(rA)
add(3, sprintf("Устойчивых признаков, набор А%s", if (A_real) "" else "*"), nA, "300–800", ok(nA >= 300 && nA <= 800))
add(8, sprintf("Полный цикл, набор А%s", if (A_real) "" else "*"), sprintf("%.1f мин; пиковая память R %.0f МБ; по модулям: %s", rA$time / 60, rA$ram,
    paste(sprintf("%d — %.0f с", rA$status$module, rA$status$time_s), collapse = "; ")), "≤ 24–48 ч", ok(rA$ok && rA$time < 24 * 3600))
if (!is.null(B)) {
  rB <- run_project("setB_full", B[1], B[2])
  nB <- if (rB$ok) nfeat(rB) else NA
  add(3, "Устойчивых признаков, набор Б", nB, "300–800", ok(!is.na(nB) && nB >= 300 && nB <= 800))
  add(8, "Полный цикл, набор Б", sprintf("%.1f мин; пиковая память R %.0f МБ", rB$time / 60, rB$ram), "≤ 24–48 ч", ok(rB$ok && rB$time < 24 * 3600))
}
stA <- rd(rA, 7, "auroc_steps.csv"); nest <- rd(rA, 7, "nested_cv.csv")
add(7, "Набор А, режим по умолчанию (Холм): признаки по этапам", paste(stA$n_features, collapse = " → "))
add(7, "Набор А, режим по умолчанию: средняя AUROC по этапам", paste(sprintf("%.3f", stA$mean_auroc_descriptive), collapse = " → "))
add(7, "Набор А, вложенная CV: точность / AUC", sprintf("%.3f / %.3f", nest$accuracy_mean, nest$auc_mean))
rC <- run_project("setA_stage1_mode", A[1], A[2], select = list(p_adjust = "stage1"),
                  filter = list(min_group_detect = 0, min_qc_detect = 0, max_qc_rsd = 1e9, blank_ratio = 0))
if (rC$ok) {
  stC <- rd(rC, 7, "auroc_steps.csv")
  add(7, "Набор А, режим совместимости с этапом 1: признаки по этапам", paste(stC$n_features, collapse = " → "), "эталон: 981 → 981 → 226 → 54 → 5",
      if (A_real) ok(abs(stC$n_features[3] - 226) / 226 < 0.1) else "эталон только для реального набора А")
  add(7, "Набор А, режим совместимости: средняя AUROC по этапам", paste(sprintf("%.3f", stC$mean_auroc_descriptive), collapse = " → "),
      "эталон: 0.541 → 0.564 → 0.642 → 0.693 → 0.873", if (A_real) "сравнить с эталоном" else "эталон только для реального набора А")
} else add(7, "Набор А, режим совместимости с этапом 1", "ошибка", "", "не пройден")

# ------------------------------------------------------------------ 4–5. набор В: пропуски и корректность
for (mp in c(0.20, 0.25, 0.30)) {
  v <- make_set_V(A[1], A[2], file.path(W, sprintf("data_setV_%02d", round(100 * mp))), n_spike = 15, missing = mp)
  r <- run_project(sprintf("setV_%02d", round(100 * mp)), v$table, v$metadata)
  if (!r$ok) { add(4, sprintf("Набор В, %d%% пропусков", round(100 * mp)), "ошибка", "", "не пройден"); next }
  im <- readRDS(file.path(dirname(v$table), "introduced_missing.rds"))
  X5 <- load_module(r$s, 5); info <- X5$info; X4 <- load_module(r$s, 4)
  # ошибка заполнения на внесённых пропусках: сравнение с истинными значениями (в шкале log2)
  tab <- data.table::fread(v$table, data.table = FALSE, check.names = FALSE, colClasses = list(character = 1))
  cols <- names(tab)[-1]; rows <- im$mask_rows
  w <- which(im$mask, arr.ind = TRUE); tv <- im$true_values
  fid <- info$id[match(cols[w[, 2]], info$original)]; smp <- load_module(r$s, 0)$meta[rows[w[, 1]], "sample"]
  okm <- !is.na(fid) & smp %in% rownames(X5$X) & !duplicated(paste(smp, fid))
  mm <- load_module(r$s, 0)$meta; single <- names(which(table(mm$sample[mm$kind == "sample"]) == 1))
  okm <- okm & smp %in% single                                        # без усреднённых повторов
  est <- X5$X[cbind(smp[okm], fid[okm])]
  hm <- half_min_impute(X4$X)[cbind(smp[okm], fid[okm])]
  nr <- function(e) sqrt(mean((log2(e) - log2(tv[okm]))^2)) / stats::sd(log2(tv[okm]))
  add(4, sprintf("Набор В, %d%% пропусков (фактически %.1f%%)", round(100 * mp), v$missing_pct),
      sprintf("обработан; NRMSE заполнения на внесённых пропусках: лес %.3f, ½ минимума %.3f (n = %d)", nr(est), nr(hm), sum(okm)),
      "без ошибок, ошибка заполнения меньше, чем у ½ минимума", ok(nr(est) < nr(hm)))
  if (mp == 0.25) {
    tr <- v$truth; tid <- info$id[match(tr$column, info$original)]
    u <- rd(r, 7, "uvf_stats.csv"); fs <- rd(r, 7, "ds_d.csv"); fr <- rd(r, 7, "ds_rfe.csv"); ne <- rd(r, 7, "nested_cv.csv")
    sel <- u$id[u$selected]
    add(5, "Набор В: найдено внесённых признаков после ОСФ", sprintf("%d из %d; ложных %d", sum(tid %in% sel), length(tid), sum(!sel %in% tid)), "≥ 80 % внесённых, ложных ≤ 1",
        ok(sum(tid %in% sel) >= 0.8 * length(tid) && sum(!sel %in% tid) <= 1))
    add(5, "Набор В: после СПФ / после РОП", sprintf("%d (внесённых %d) / %d (внесённых %d)", ncol(fs) - 2, sum(tid %in% names(fs)), ncol(fr) - 2, sum(tid %in% names(fr))))
    add(5, "Набор В: вложенная CV, точность / AUC", sprintf("%.3f / %.3f", ne$accuracy_mean, ne$auc_mean))
    data.table::fwrite(data.frame(tr, id = tid, found_uvf = tid %in% sel, in_ds_d = tid %in% names(fs), in_ds_rfe = tid %in% names(fr)),
                       file.path(OUT, "setV_truth_vs_found.csv"), bom = TRUE)
  }
}

# ------------------------------------------------------------------ 6. воспроизводимость
keyf <- c("05_impute/imputed.csv", "06_correct/corrected_log2.csv", "07_select/uvf_stats.csv", "07_select/ds.csv", "07_select/ds_d.csv",
          "07_select/ds_rfe.csv", "07_select/model_accuracy.csv", "07_select/nested_cv.csv")
h <- lapply(1:3, function(i) { r <- run_project(sprintf("repro_same_seed_%d", i), A[1], A[2], project = list(cores = if (i == 3) max(2, CORES) else CORES)); md5(r, keyf) })
same <- all(vapply(h[-1], function(x) identical(x, h[[1]]), TRUE))
add(6, "3 запуска с одним seed (третий — в параллельном режиме)", if (same) "все ключевые таблицы побайтно совпадают (md5)" else "есть различия",
    "полное совпадение", ok(same))
sets <- lapply(1:5, function(i) { r <- run_project(sprintf("repro_seed_%d", i), A[1], A[2], project = list(seed = 1000 + i, cores = CORES))
  list(rfe = setdiff(names(rd(r, 7, "ds_rfe.csv")), c("sample", "Label")), uvf = setdiff(names(rd(r, 7, "ds.csv")), c("sample", "Label"))) })
stab <- function(k) { L <- lapply(sets, `[[`, k); length(Reduce(intersect, L)) / max(1, length(Reduce(union, L))) }
add(6, "5 запусков с разными seed: доля признаков, отобранных во всех запусках (РОП)", sprintf("%.2f (во всех %d, хотя бы в одном %d)", stab("rfe"),
    length(Reduce(intersect, lapply(sets, `[[`, "rfe"))), length(Reduce(union, lapply(sets, `[[`, "rfe")))))
add(6, "5 запусков с разными seed: то же для ОСФ (не зависит от seed)", sprintf("%.2f", stab("uvf")), "1.00", ok(stab("uvf") == 1))

# ------------------------------------------------------------------ 9. полнота отчёта
rep <- file.path(rA$dir, "09_report", "report.html"); html <- if (file.exists(rep)) paste(readLines(rep, encoding = "UTF-8"), collapse = "\n") else ""
need <- c("Сведения о прогоне", "0. Входные данные", "2. Контроль качества", "3. Фильтр", "5. Заполнение", "6. Коррекция", "7. Статистика",
          "8. Визуализация", "9. База признаков", "Журнал прогона", "sessionInfo")
miss <- need[!vapply(need, grepl, TRUE, x = html, fixed = TRUE)]
nimg <- lengths(regmatches(html, gregexpr("data:image/png;base64", html)))
add(9, "HTML-отчёт: разделы и рисунки", sprintf("разделов %d из %d; встроенных рисунков %d; размер %.1f МБ", length(need) - length(miss), length(need), nimg,
    file.size(rep) / 2^20), "все разделы, все рисунки", ok(!length(miss) && nimg >= 8))
add(9, "Табличные результаты", paste(length(list.files(rA$dir, pattern = "\\.csv$", recursive = TRUE)), "CSV-файлов; база:", paste(list.files(file.path(rA$dir, "features_db")), collapse = ", ")))
file.copy(rep, file.path(OUT, "report_setA.html"), overwrite = TRUE)

# ------------------------------------------------------------------ 10. установка с нуля
add(10, "Установка на чистый ПК с Windows", "выполняется вручную по чек-листу (docs/user_guide.md, раздел «Установка»)", "время установки, проблемы", "не выполнен в среде разработки")

writeLines(c(paste("UroMetabR", UROMETABR_VERSION), utils::capture.output(utils::sessionInfo())), file.path(OUT, "session_info_ispytaniya.txt"))
cat("\nПротокол:", file.path(OUT, "protocol.csv"), "\n")

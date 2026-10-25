# Функциональные испытания на контрольных наборах с известной истиной (tests/make_control_data.R).

root <- getOption("urometabr.root")
ctl_cfg <- function(dir, out, perm = 0) {
  cfg <- read_config(file.path(root, "config", "config_example.yaml"), quick = TRUE)
  cfg$input$data <- file.path(dir, "feature_table.csv"); cfg$input$metadata <- file.path(dir, "metadata.csv")
  cfg$output <- out; cfg$multivariate$n_perm <- perm; cfg$report$fig_dpi <- 72
  cfg
}
series <- file.path(root, "example", "contaminant_series.csv")
d1 <- file.path(tempdir(), "ctl1"); unlink(d1, recursive = TRUE)
make_control_data(d1, seed = 11, series_file = series)
truth <- read.csv(file.path(d1, "truth.csv"), stringsAsFactors = FALSE); rownames(truth) <- truth$feature
ts <- read.csv(file.path(d1, "truth_samples.csv"), colClasses = c(sample = "character"))
invisible(capture.output(R1 <- run_pipeline(ctl_cfg(d1, file.path(tempdir(), "out1")))))

test_that("фильтры удаляют нестабильные, редкие и фоновые признаки", {
  kept <- rownames(R1$features)
  tt <- truth[kept, "type"]
  expect_lte(mean(truth$feature[truth$type == "noisy"] %in% kept), 0.15)   # ≥ 85 % шумных удалено
  expect_equal(sum(tt == "rare"), 0)
  expect_equal(sum(tt == "blank"), 0)
  expect_gte(mean(truth$feature[truth$type %in% c("BC_up", "BC_down")] %in% kept), 0.9)
})

test_that("PQN восстанавливает истинные факторы разведения", {
  f <- R1$meta[ts$sample, "dilution"]
  expect_gt(cor(log(f), log(ts$dilution_true)), 0.97)
})

test_that("аддукты объединяются с исходными ионами", {
  ad <- truth[truth$type == "adduct" & truth$feature %in% rownames(R1$features) & truth$parent %in% rownames(R1$features), ]
  same <- R1$features[ad$feature, "compound"] == R1$features[ad$parent, "compound"]
  s <- R1$meta$sample[R1$meta$study]
  X <- data.table::fread(file.path(d1, "feature_table.csv"), data.table = FALSE, check.names = FALSE)
  miss <- colMeans(X[X$Sample %in% s, ad$feature, drop = FALSE] == 0)
  expect_gt(mean(same[miss <= 0.05]), 0.95)       # аддукты без пропусков объединяются почти всегда
  expect_gt(mean(same), 0.7)                       # с пропусками (импутация) корреляция ниже
  message(sprintf("аддукты: объединено %.0f%% (без пропусков %.0f%%, n = %d)", 100 * mean(same), 100 * mean(same[miss <= 0.05]), length(same)))
})

test_that("индекс загрязнения восстанавливает истинную экспозицию", {
  expect_gt(cor(R1$meta[ts$sample, "contam_index"], ts$contam_true, method = "spearman"), 0.9)
})

test_that("находки: истинные эффекты найдены, ложных мало, экспозиции и возраст исключены", {
  r <- R1$res$BC_vs_Control
  f <- rownames(r)[r$finding]; tt <- truth[f, "type"]
  expect_gte(sum(tt %in% c("BC_up", "BC_down")), 20)                   # чувствительность ≥ 40 % при эффекте 0.7 log2
  expect_lte(mean(!tt %in% c("BC_up", "BC_down", "drug_BC")), 0.10)    # доля ложных находок ≤ 10 %
  expect_true(all(sign(r[f[tt == "BC_up"], "beta_model"]) > 0))
  expect_false(any(tt %in% c("age_only", "series")))
  expect_false(any(f %in% truth$feature[truth$type == "drug_BC" & !grepl("low score", truth$feature)]))
  expect_equal(sum(R1$res$PC_vs_Control$finding), 0)                   # у PC истинных эффектов нет
})

test_that("поправка на возраст устраняет ложные отличия, вызванные возрастом", {
  r <- R1$res$PC_vs_Control; age <- intersect(rownames(r), truth$feature[truth$type == "age_only"])
  expect_gt(sum(r[age, "q_MW"] < 0.05), 0)            # без поправки возрастные признаки «значимы»
  expect_equal(sum(r[age, "q_model"] < 0.05), 0)      # в модели с возрастом — нет
})

test_that("набор без групповых эффектов: нет находок", {
  d0 <- file.path(tempdir(), "ctl0"); unlink(d0, recursive = TRUE)
  make_control_data(d0, seed = 12, null = TRUE, series_file = series)
  cfg <- ctl_cfg(d0, file.path(tempdir(), "out0")); cfg$multivariate$run <- FALSE
  invisible(capture.output(R0 <- run_pipeline(cfg)))
  expect_equal(sum(vapply(R0$res, function(r) sum(r$finding), 0)), 0)
})

test_that("кросс-валидация без утечки: случайные метки дают AUC ≈ 0.5", {
  cfg <- ctl_cfg(d1, tempdir(), perm = 50); cfg$multivariate$cv_repeats <- 5
  ctrl <- R1$meta$sample[R1$meta$study & R1$meta$group == "Control"]
  set.seed(99); a <- sample(ctrl, 22); b <- setdiff(ctrl, a)
  e <- evaluate_models(R1$Lc, R1$meta, a, b, colnames(R1$Lc), "random", 1, cfg)
  pnl <- e$table[grepl("^панель top-5 \\+", e$table$model), ]
  expect_lt(abs(pnl$AUC_mean - 0.5), 0.15)
  expect_gt(pnl$perm_p, 0.05)
  # контроль: отбор признаков ВНЕ CV (ошибка макета 1 этапа) даёт завышенную AUC на тех же случайных метках
  y <- c(rep(1, length(a)), rep(0, length(b))); X <- R1$Lc[c(a, b), ]
  sel <- order(mw_p_fast(X, y))[1:5]
  set.seed(1); cv <- cv_predict(X[, sel], y, "covariate", 5, 5)
  expect_gt(auc_score(y, cv$pred), pnl$AUC_mean + 0.15)
})

test_that("воспроизводимость: повторный запуск и параллельный режим дают идентичные таблицы", {
  cfg <- ctl_cfg(d1, file.path(tempdir(), "rep_a"), perm = 20); cfg$cores <- 1
  invisible(capture.output(run_pipeline(cfg)))
  cfg$output <- file.path(tempdir(), "rep_b"); cfg$cores <- 2
  invisible(capture.output(run_pipeline(cfg)))
  a <- read.csv(file.path(tempdir(), "rep_a", "manifest_md5.csv")); b <- read.csv(file.path(tempdir(), "rep_b", "manifest_md5.csv"))
  a <- a[a$file != "log.csv", ]; b <- b[b$file != "log.csv", ]       # журнал содержит время обработки
  expect_equal(a, b)
  expect_gte(nrow(a), 10)
})

test_that("полнота результатов: все таблицы, рисунки и отчёт сформированы", {
  out <- file.path(tempdir(), "out1")
  for (f in c("results.xlsx", "report.html", "session_info.txt", "manifest_md5.csv")) expect_true(file.exists(file.path(out, f)))
  for (t in c("summary_univariate", "findings", "multivariate_CV", "features", "samples", "demographics", "presence_absence",
              "normalization_sensitivity", "age_matched_pairs", "processed_log2_matrix"))
    expect_true(file.exists(file.path(out, "tables", paste0(t, ".csv"))), info = t)
  figs <- list.files(file.path(out, "figures"), pattern = "png$")
  expect_gte(length(figs), 19)
  expect_true(all(file.size(file.path(out, "figures", figs)) > 5000))
  expect_gt(file.size(file.path(out, "report.html")), 1e5)
})

# Модульные тесты: каждая функция сверяется с эталонной реализацией из базового R / CRAN.

test_that("Манн–Уитни совпадает с wilcox.test (связи, поправка на непрерывность)", {
  set.seed(1)
  A <- matrix(round(rnorm(200, 1), 1), 20); B <- matrix(round(rnorm(150), 1), 15)
  r <- mw_matrix(A, B)
  ref <- vapply(1:10, function(j) suppressWarnings(wilcox.test(A[, j], B[, j], exact = FALSE, correct = TRUE)$p.value), 0)
  expect_equal(r$p, ref, tolerance = 1e-10)
  expect_equal(r$auc, vapply(1:10, function(j) as.numeric(pROC::auc(pROC::roc(c(rep(1, 20), rep(0, 15)), c(A[, j], B[, j]),
                                                                            direction = "<", quiet = TRUE))), 0), tolerance = 1e-10)
})

test_that("модель с ковариатами совпадает с lm()", {
  set.seed(2)
  n <- 40; meta <- data.frame(sample = paste0("s", 1:n), sex = sample(c("M", "F"), n, TRUE), age = runif(n, 30, 80),
                              contam_index = rnorm(n), row.names = paste0("s", 1:n))
  meta$age[3] <- NA
  L <- matrix(rnorm(n * 5), n, dimnames = list(meta$sample, paste0("f", 1:5)))
  case <- meta$sample[1:18]; ref <- meta$sample[19:40]
  r <- ols_group(L, meta, case, ref, c("sex", "age", "contam_index"))
  for (j in 1:5) {
    d <- data.frame(y = L[, j], g = as.numeric(meta$sample %in% case), sex = as.numeric(meta$sex == "M"), age = meta$age, ci = meta$contam_index)
    s <- summary(lm(y ~ g + sex + age + ci, d))$coefficients["g", ]
    expect_equal(unname(r$beta[j]), unname(s["Estimate"]), tolerance = 1e-8)
    expect_equal(unname(r$p[j]), unname(s["Pr(>|t|)"]), tolerance = 1e-8)
  }
  expect_equal(r$n_case + r$n_ref, 39)          # complete case: одна проба без возраста исключена
})

test_that("BH, AUC и логистическая регрессия совпадают с эталонами", {
  p <- c(0.001, 0.02, NA, 0.04, 0.5); expect_equal(bh(p), p.adjust(p, "BH"))
  set.seed(3); y <- rep(0:1, 30); s <- rnorm(60) + y
  expect_equal(auc_score(y, s), as.numeric(pROC::auc(pROC::roc(y, s, direction = "<", quiet = TRUE))))
  X <- matrix(rnorm(120), 60)
  w <- fit_logreg(X, y, C = 1e8)
  expect_equal(w, unname(coef(glm(y ~ X, family = binomial))), tolerance = 1e-5)
})

test_that("PQN восстанавливает известные факторы разведения, MSTUS — тоже", {
  set.seed(4)
  base <- 2^rnorm(300, 15, 2); f <- exp(rnorm(30, 0, 0.6))
  X <- outer(f, base) * 2^matrix(rnorm(9000, 0, 0.1), 30); rownames(X) <- c(paste0("QC", 1:5), paste0("s", 1:25))
  X[sample(length(X), 300)] <- NA
  meta <- data.frame(kind = rep(c("QC", "sample"), c(5, 25)), row.names = rownames(X))
  pq <- normalize(X, meta, "PQN")$f; ms <- normalize(X, meta, "MSTUS")$f
  expect_gt(cor(log(pq), log(f)), 0.99)
  expect_gt(cor(log(ms), log(f)), 0.99)
  expect_equal(median(pq[1:5] / f[1:5] * median(f[1:5]) ), median(f[1:5]) * median(pq[1:5] / f[1:5]))
})

test_that("импутация ½ минимума и log2", {
  X <- matrix(c(4, NA, 16, 8, 2, NA), 3); r <- impute_log(X)
  expect_equal(r[2, 1], log2(2)); expect_equal(r[3, 2], log2(1)); expect_false(anyNA(r))
})

test_that("фильтры детекции, бланков и RSD", {
  cfg <- default_config(); cfg$design$groups <- c("A", "B")
  X <- cbind(good = c(10, 11, 10, 10, 9, 11, NA, NA), noisy = c(1, 20, 5, 10, 9, 11, NA, NA),
             blank = c(10, 10, 10, 10, 10, 10, 9, 9), rare = c(10, 10, NA, NA, NA, NA, NA, 10))
  meta <- data.frame(kind = c("QC", "QC", "sample", "sample", "sample", "sample", "blank", "blank"),
                     group = c("QC", "QC", "A", "A", "B", "B", "Blank", "Blank"), study = c(F, F, T, T, T, T, F, F))
  invisible(capture.output(r <- filter_features(X, meta, cfg)))
  expect_equal(colnames(r$X), "good")
  expect_false(r$table["blank", "pass_detection"]); expect_false(r$table["rare", "pass_detection"])
  expect_true(r$table["noisy", "RSD_QC_raw"] > 0.3)
})

test_that("подбор пар: 1:1, внутри пола, в пределах допуска, оптимально", {
  meta <- data.frame(sample = c(paste0("c", 1:4), paste0("r", 1:5)), sex = c("M", "M", "F", "M", "M", "M", "F", "F", "M"),
                     age = c(60, 70, 50, 90, 61, 69, 52, 30, 75))
  rownames(meta) <- meta$sample
  cfg <- default_config()
  p <- age_match(meta, meta$sample[1:4], meta$sample[5:9], cfg)
  expect_equal(nrow(p), 3)                                        # c4 (90 лет) без пары
  expect_true(all(p$diff <= 5)); expect_false(any(duplicated(p$ref)))
  expect_true(all(meta[p$case, "sex"] == meta[p$ref, "sex"]))
  expect_equal(sum(p$diff), 1 + 1 + 2)                            # оптимальное назначение
})

test_that("масса формулы и группировка аддуктов", {
  expect_equal(formula_mass("C6H12O6"), 180.0633881, tolerance = 1e-6)
  set.seed(5); n <- 30; a <- rnorm(n, 15); 
  L <- cbind(f1 = a, f2 = a + rnorm(n, 0, 0.05), f3 = rnorm(n, 15), f4 = rnorm(n, 15))
  info <- data.frame(feature = colnames(L), rt_min = c(5, 5.01, 5.005, 9), mz = c(200, 200 + 21.981942, 300, 200), row.names = colnames(L))
  g <- group_features(L, info, stats::setNames(c(10, 3, 5, 5), colnames(L)))
  expect_equal(g["f1", "compound"], g["f2", "compound"]); expect_false(g["f3", "compound"] == g["f1", "compound"])
  expect_true(g["f1", "is_representative"]); expect_match(g["f2", "ion_relations"], "Na-H")
})

test_that("чтение экспорта MS-DIAL даёт ту же матрицу, что и CSV", {
  d <- file.path(tempdir(), "msd"); unlink(d, recursive = TRUE)
  make_control_data(d, seed = 7, n_feat = 600, series_file = file.path(getOption("urometabr.root"), "example", "contaminant_series.csv"))
  write_msdial_export(d)
  cfg <- default_config(); cfg$input$data <- file.path(d, "feature_table.csv")
  w <- read_wide(cfg$input$data, cfg); m <- read_msdial(file.path(d, "msdial_area.txt"), cfg)
  expect_equal(unname(m$X[rownames(w$X), ]), unname(w$X), tolerance = 1e-6)
  expect_equal(unname(m$kind[rownames(w$X)]), unname(w$kind))
  expect_equal(m$info$mz, w$info$mz)
})

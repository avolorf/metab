# Модульные тесты алгоритмов и исправлений кода этапа 1.
root <- getOption("urometabr.root")
code <- unlist(lapply(c(list.files(file.path(root, "R"), full.names = TRUE), file.path(root, "app", "app.R"), file.path(root, "run_analysis.R")),
                      readLines, encoding = "UTF-8"))
code <- code[!grepl("^\\s*#", code)]

test_that("исправление 1: нет жёстко прописанных путей и setwd()", {
  expect_false(any(grepl("setwd\\(", code)))
  expect_false(any(grepl("[\"'][A-Za-z]:/", code)))
})

test_that("исправление 5: метки групп — из файла метаданных, а не из имени файла", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("file,group,injection_order,replicate,sample", "a 1 x,BC,2,1,A", "b_2,CG,1,1,B", "a 1 x_r2,BC,3,2,A", "QC1,QC,4,1,QC1"), f)
  m <- read_metadata(f)
  expect_equal(unname(m["a 1 x", "group"]), "BC"); expect_equal(m$kind, c("sample", "sample", "sample", "QC"))
  expect_equal(m["a 1 x_r2", "replicate"], 2L)
  writeLines(c("file,grp", "a,BC"), f); expect_error(read_metadata(f), "обязательных столбцов")
})

test_that("исправление 6: поправка — параметр; по умолчанию Холм, применяется по всем признакам", {
  expect_equal(default_settings()$select$p_adjust, "holm")
  set.seed(1); L <- matrix(rnorm(40 * 200), 40); colnames(L) <- sprintf("F%03d", 1:200); g <- rep(c("A", "B"), 20)
  L[g == "B", 1:5] <- L[g == "B", 1:5] + 2
  h <- uvf(L, g, "holm"); b <- uvf(L, g, "BH"); n <- uvf(L, g, "none")
  expect_equal(h$p_adj, p.adjust(h$p, "holm")); expect_equal(b$p_adj, p.adjust(b$p, "BH"))
  expect_true(all(h$id[h$selected] %in% sprintf("F%03d", 1:5)))
  expect_gt(sum(n$selected), sum(h$selected))                       # без поправки отбираются и случайные признаки
  ref <- t.test(L[g == "A", 1], L[g == "B", 1], var.equal = h$test[1] == "Стьюдент")$p.value
  if (h$test[1] != "Уилкоксон") expect_equal(h$p[1], ref)
})

test_that("ОСФ: выбор теста (Уилкоксон / Уэлч / Стьюдент) и режим совместимости с этапом 1", {
  set.seed(2); g <- rep(c("A", "B"), each = 30)
  L <- cbind(F1 = c(rexp(30), rexp(30) + 1), F2 = c(rnorm(30, 0, 1), rnorm(30, 1, 5)), F3 = c(rnorm(30), rnorm(30, 1)))
  u <- uvf(L, g, "holm")
  expect_equal(u$test, c("Уилкоксон", "Уэлч", "Стьюдент"))
  expect_equal(u$p[1], wilcox.test(L[1:30, 1], L[31:60, 1], exact = FALSE)$p.value)
  s1 <- uvf(L, g, "stage1"); expect_equal(s1$p_adj, s1$p)
})

test_that("исправление 7: частота встречаемости считается через table()", {
  expect_false(any(grepl("tuplicated\\(|library\\(tuple\\)", code)))
  expect_true(any(grepl("table\\(unlist\\(tops\\)\\)", code)))
})

test_that("исправление 3: ds_rfe собирается из predictors(rfe)", {
  set.seed(3); g <- factor(rep(c("A", "B"), 25)); X <- matrix(rnorm(50 * 6), 50); colnames(X) <- paste0("F", 1:6)
  X[g == "B", 1:2] <- X[g == "B", 1:2] + 2
  s <- default_settings(); s$select$rfe_folds <- 3; s$select$rfe_repeats <- 1
  r <- run_rfe(data.frame(Label = g, X), s, 1)
  expect_true(all(c("F1", "F2") %in% caret::predictors(r)))
  expect_true(any(grepl("f_rfe <- caret::predictors\\(out\\$rfe\\)", code)))
})

test_that("исправление 8: цвета групп берутся из данных (нет ручных ключей с кириллицей и пробелами)", {
  col <- group_colors(c("CG", "BC", "QC")); expect_equal(names(col), c("BC", "CG", "QC"))
  expect_true(all(grepl("^#[0-9A-F]{6}$", col)))
  expect_false(any(grepl("\" #[0-9a-fA-F]{6}\"", code)))
})

test_that("исправления 9–10: один режим параллельных вычислений, minFraction — параметр", {
  x <- readLines(file.path(root, "R", "m1_peaktable.R"), encoding = "UTF-8")
  expect_equal(sum(grepl("BiocParallel::register\\(", x)), 1)
  expect_true(any(grepl("minFraction = x\\$min_fraction", x)))
  expect_equal(default_settings()$xcms$min_fraction, 0.5)
})

test_that("EigenMS: удаляет внесённый тренд и сохраняет различия групп", {
  set.seed(4); n <- 40; p <- 300; g <- rep(c("A", "B"), each = 20)
  L <- matrix(rnorm(n * p, 15, 0.3), n); L[g == "B", 1:20] <- L[g == "B", 1:20] + 1
  trend <- sin(seq(0, 3 * pi, length.out = n)); load <- rnorm(p, 0, 0.8)
  Lb <- L + outer(trend, load)
  e <- eigenms(Lb, g, n_perm = 100, seed = 1)
  expect_gte(e$n_trends, 1)
  expect_lt(abs(cor(e$L[, 30] - L[, 30], trend)), abs(cor(Lb[, 30] - L[, 30], trend)))
  expect_lt(mean(abs(e$L - L - (colMeans(e$L - L)[col(L)]))), mean(abs(Lb - L - colMeans(Lb - L)[col(L)])) / 3)
  gm <- function(M) colMeans(M[g == "B", ]) - colMeans(M[g == "A", ])
  expect_equal(gm(e$L), gm(Lb), tolerance = 1e-8)                  # различия групп не меняются
  e0 <- eigenms(L, g, n_perm = 100, seed = 1); expect_equal(e0$n_trends, 0)
})

test_that("заполнение пропусков случайным лесом: ошибка меньше, чем у ½ минимума и среднего", {
  set.seed(5); n <- 80; z <- matrix(rnorm(n * 3), n)
  L <- cbind(z %*% matrix(rnorm(3 * 30), 3), matrix(rnorm(n * 10), n)) + 15
  X <- 2^L; colnames(X) <- paste0("F", 1:40); M <- X; idx <- sample(length(M), 0.15 * length(M)); M[idx] <- NA
  r <- rf_impute(M, maxiter = 5, ntree = 100, seed = 1)
  err <- function(Y) sqrt(mean((log2(Y[idx]) - L[idx])^2))
  mean_imp <- apply(M, 2, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v })
  expect_false(anyNA(r$X)); expect_lt(err(r$X), err(half_min_impute(M))); expect_lt(err(r$X), err(mean_imp))
  r2 <- rf_impute(M, maxiter = 5, ntree = 100, seed = 1); expect_identical(r$X, r2$X)   # исправление 11: воспроизводимость
})

test_that("усреднение повторов и сортировка по порядку ввода (модуль 4)", {
  d <- file.path(tempdir(), "m4"); unlink(d, recursive = TRUE); dir.create(d)
  s <- default_settings(); s$project$dir <- d
  X <- matrix(c(1, 3, 10, 5, NA, 7, 2, 4, 20, 6, 8, 9), 6, dimnames = list(c("a1", "a2", "b", "c", "QC1", "QC2"), c("F1", "F2")))
  meta <- data.frame(file = rownames(X), sample = c("A", "A", "B", "C", "QC1", "QC2"), group = c("BC", "BC", "CG", "CG", "QC", "QC"),
                     kind = c(rep("sample", 4), "QC", "QC"), injection_order = c(5, 6, 1, 3, 2, 4), row.names = rownames(X))
  save_module(s, 3, list(X = X, info = data.frame(id = c("F1", "F2")), meta = meta))
  invisible(capture.output(module4_samples(s))); r <- load_module(s, 4)
  expect_equal(rownames(r$X), c("B", "C", "A")); expect_equal(unname(r$X["A", "F1"]), 2); expect_equal(unname(r$X["A", "F2"]), 3)
})

test_that("модуль 0: понятные ошибки на неверных входных данных", {
  d <- file.path(tempdir(), "m0"); unlink(d, recursive = TRUE); dir.create(d)
  writeLines(c("file,F1,F2", "x1,1,2", "x2,3,4"), file.path(d, "t.csv")); writeLines(c("file,group", "y1,A", "y2,B"), file.path(d, "m.csv"))
  s <- read_settings(new_project(d, "t.csv", "m.csv"))
  expect_error(module0_input(s), "Слишком мало совпадений")
  s$input$table_file <- "нет.csv"; expect_error(module0_input(s), "Не найден файл таблицы")
  s$input$type <- "raw"; s$input$raw_dir <- "нет_папки"; expect_error(module0_input(s), "Не найдена папка")
  write_raw_stubs(file.path(d, "raw"), c("y1", "y2")); writeLines("garbage", file.path(d, "raw", "bad.mzML"))
  s$input$raw_dir <- "raw"; expect_error(module0_input(s), "bad.mzML")
})

test_that("модуль 3: стандартный и строгий профили фильтра", {
  d <- file.path(tempdir(), "m3"); unlink(d, recursive = TRUE); dir.create(d)
  s <- default_settings(); s$project$dir <- d
  set.seed(6); n <- 12; X <- matrix(2^rnorm(n * 4, 15, 0.05), n, dimnames = list(c(paste0("s", 1:8), paste0("QC", 1:4)), paste0("F", 1:4)))
  X[1:4, "F2"] <- NA                                  # только в одной группе — стандартный да, строгий нет
  X[9:12, "F3"] <- X[9:12, "F3"] * c(1, 1.3, 0.8, 1.2) # RSD в QC ≈ 20 % — стандартный да, строгий нет
  X[1:7, "F4"] <- NA                                  # почти нигде — оба нет
  meta <- data.frame(file = rownames(X), group = c(rep(c("A", "B"), each = 4), rep("QC", 4)), kind = c(rep("sample", 8), rep("QC", 4)),
                     row.names = rownames(X))
  save_module(s, 2, list(X = X, info = data.frame(id = colnames(X), mz = 1:4, rt_min = 1:4), meta = meta))
  invisible(capture.output(module3_filter(s))); r <- load_module(s, 3)
  expect_equal(r$standard_ids, c("F1", "F2", "F3")); expect_equal(r$strict_ids, "F1")
  expect_true(all(r$strict_ids %in% r$standard_ids)); expect_equal(colnames(r$X), r$standard_ids)
  expect_setequal(unique(r$steps$profile), c("стандартный", "строгий"))
  s$filter$analysis_profile <- "strict"; s$filter$strict$max_qc_rsd <- 0.3
  invisible(capture.output(module3_filter(s))); expect_equal(colnames(load_module(s, 3)$X), c("F1", "F3"))
})

test_that("модуль 9: пустой строгий профиль не останавливает обработку", {
  d <- file.path(tempdir(), "m9"); unlink(d, recursive = TRUE); dir.create(d)
  f <- new_project(d); s <- read_settings(f)
  L <- matrix(rnorm(40), 10, dimnames = list(paste0("s", 1:10), paste0("F", 1:4)))
  info <- data.frame(id = colnames(L), mz = 1:4, rt_min = 1:4, name = "", original = colnames(L))
  meta <- data.frame(sample = rownames(L), group = rep(c("A", "B"), 5), kind = "sample", row.names = rownames(L))
  save_module(s, 3, list(X = 2^L, info = info, meta = meta, filter = data.frame(id = colnames(L), rsd_qc = 0.1),
                         strict_ids = character(0), standard_ids = colnames(L)))
  save_module(s, 4, list(X = 2^L, info = info, meta = meta))
  save_module(s, 7, list(L = L, L_raw = L, info = info, meta = meta,
                         sel = list(uvf = data.frame(id = colnames(L), p = 0.5, p_adj = 1), f_uvf = character(0), f_spf = character(0), f_rfe = character(0))))
  log_start(s); s$report$rmarkdown <- FALSE
  expect_error(invisible(capture.output(module9_report(s))), NA)
  db <- read.csv(file.path(d, "features_db", "features_db.csv"))
  expect_equal(nrow(db), 4); expect_true(all(grepl("strict→standard", db$db_profile)))
})

test_that("отбор по умолчанию — на данных без EigenMS", {
  expect_equal(default_settings()$select$selection_input, "imputed")
  expect_true(any(grepl('Ls <- if \\(sl\\$selection_input == "imputed"\\) d\\$L_raw else d\\$L', code)))
})

test_that("чтение таблицы: интенсивности больше 2^31 не искажаются (integer64)", {
  f <- tempfile(fileext = ".csv"); md <- data.frame(file = c("s1", "s2"), group = c("BC", "CG"))
  writeLines(c("file,F1,F2", "s1,4212240040,10", "s2,892085242,20"), f)
  X <- read_feature_table(f, md)$X
  expect_equal(unname(X[, "F1"]), c(4212240040, 892085242))
})

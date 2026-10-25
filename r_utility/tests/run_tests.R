#!/usr/bin/env Rscript
# Проверка корректности UroMetabR: модульные тесты алгоритмов и функциональные испытания
# на контрольных наборах данных с известной истиной.
#   Rscript tests/run_tests.R            (из каталога r_utility)
if (!isTRUE(l10n_info()$`UTF-8`)) invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
root <- normalizePath(".")
options(urometabr.dir = file.path(root, "R"), urometabr.root = root)
for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
source(file.path(root, "run_analysis.R"), encoding = "UTF-8")
source(file.path(root, "tests", "make_control_data.R"), encoding = "UTF-8")
res <- testthat::test_dir(file.path(root, "tests", "testthat"), reporter = testthat::ProgressReporter$new(show_praise = FALSE),
                          stop_on_failure = FALSE, load_package = "none", env = globalenv())
df <- as.data.frame(res)
out <- data.frame(file = df$file, test = df$test, n_checks = df$nb, failed = df$failed, error = df$error, time_s = round(df$real, 1))
data.table::fwrite(out, file.path(root, "tests", "test_results.csv"))
cat(sprintf("\nИтого: тестов %d, проверок %d, не пройдено %d, ошибок %d\n", nrow(out), sum(out$n_checks), sum(out$failed), sum(out$error)))
quit(status = as.integer(sum(out$failed) + sum(out$error) > 0))

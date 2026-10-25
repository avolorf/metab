#!/usr/bin/env Rscript
# Модульные тесты: Rscript tests/run_tests.R (из папки утилиты)
if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows")
  for (loc in c("C.UTF-8", "en_US.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
root <- normalizePath(".")
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8"); load_urometabr(root)
source(file.path(root, "tests", "make_test_data.R"), encoding = "UTF-8")
res <- as.data.frame(testthat::test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = FALSE,
                                        load_package = "none", env = globalenv()))
out <- data.frame(test = res$test, checks = res$nb, failed = res$failed, error = res$error)
data.table::fwrite(out, file.path(root, "tests", "unit_test_results.csv"), bom = TRUE)
cat(sprintf("\nТестов %d, проверок %d, не пройдено %d, ошибок %d\n", nrow(out), sum(out$checks), sum(out$failed), sum(out$error)))
quit(status = as.integer(sum(out$failed) + sum(out$error) > 0))

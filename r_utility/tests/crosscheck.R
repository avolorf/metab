#!/usr/bin/env Rscript
# Сверка UroMetabR с эталонной реализацией алгоритма на Python (pipeline/ репозитория) на контрольном наборе.
#   Rscript tests/crosscheck.R   (из каталога r_utility; нужен python3 с numpy, pandas, scipy, statsmodels)
if (!isTRUE(l10n_info()$`UTF-8`)) invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
options(urometabr.dir = "R")
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
source("tests/make_control_data.R", encoding = "UTF-8")
d <- file.path(tempdir(), "xc"); out <- file.path(tempdir(), "xc_py")
make_control_data(d, seed = 21, series_file = "example/contaminant_series.csv")
stopifnot(system2("python3", c("tests/crosscheck_python.py", d, out)) == 0)
cfg <- read_config("config/config_example.yaml"); cfg$input$data <- file.path(d, "feature_table.csv"); cfg$input$metadata <- file.path(d, "metadata.csv")
invisible(capture.output({
  ft <- read_feature_table(cfg); meta <- load_metadata(cfg$input$metadata, ft$kind, cfg)
  study <- meta$sample[meta$study]; use <- meta$sample[meta$study | meta$kind %in% c("QC", "blank")]
  flt <- filter_features(ft$X[use, ], meta[use, ], cfg); X1 <- flt$X
  ana <- meta$sample[meta$study | meta$kind == "QC"]
  n <- normalize(X1[ana, ], meta, "PQN"); L <- impute_log(n$X); fm <- normalize(X1[ana, ], meta, "MSTUS")$f
  g <- group_features(L[study, ], ft$info[colnames(L), ], col_medians(X1[study, ]))
  ser <- match_series(ft$info[colnames(L), ], cfg$input$contaminant_series)
  ci <- contamination_index(L, ser, study); meta$contam_index <- NA; meta[names(ci$index), "contam_index"] <- ci$index
  cm <- list(name = "BC_vs_Control", case = "BC", ref = "Control")
  s <- comparison_sets(meta, cm); Ls <- L[study, ]
  r1 <- mw(Ls, s$case, s$ref); r2 <- ols_group(Ls, meta, s$case, s$ref, c("sex", "age", "contam_index"))
  pr <- age_match(meta, s$case, s$ref, cfg)
}))
pyf <- read.csv(file.path(out, "py_filter.csv"), row.names = 1, check.names = FALSE)
pyk <- rownames(pyf)[pyf$kept == "True"]
fac <- read.csv(file.path(out, "py_factors.csv"), row.names = 1); rownames(fac) <- as.character(rownames(fac))
st <- read.csv(file.path(out, "py_stats.csv"), row.names = 1, check.names = FALSE)
gr <- read.csv(file.path(out, "py_groups.csv"), row.names = 1, check.names = FALSE)
pc <- read.csv(file.path(out, "py_contam.csv"), row.names = 1)
pp <- read.csv(file.path(out, "py_pairs.csv"))
same_partition <- function(a, b) { t <- table(a, b); all(rowSums(t > 0) == 1) && all(colSums(t > 0) == 1) }
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
res <- data.frame(
  check = c("набор признаков после фильтров", "RSD в QC (макс. отн. расхождение)", "факторы PQN", "факторы MSTUS",
            "разбиение признаков на соединения", "индекс загрязнения", "p Манна–Уитни", "AUC", "log2FC модели", "p модели",
            "число пар по возрасту", "сумма |Δвозраст| в парах"),
  value = c(sprintf("%d / %d, совпадение %s", length(colnames(X1)), length(pyk), setequal(colnames(X1), pyk)),
            sprintf("%.2e", rel(flt$table[pyk, "RSD_QC_raw"], pyf[pyk, "RSD_QC_raw"])),
            sprintf("%.2e", rel(n$f[rownames(fac)], fac$pqn)), sprintf("%.2e", rel(fm[rownames(fac)], fac$mstus)),
            as.character(same_partition(g[rownames(gr), "compound"], gr$compound)),
            sprintf("%.2e (проб: %d)", max(abs(meta[rownames(pc), "contam_index"] - pc$contam_index), na.rm = TRUE), sum(!is.na(pc$contam_index))),
            sprintf("%.2e", rel(r1$p, st[colnames(Ls), "p_MW"])), sprintf("%.2e", max(abs(r1$auc - st[colnames(Ls), "AUC"]))),
            sprintf("%.2e", max(abs(r2$beta - st[colnames(Ls), "beta"]))), sprintf("%.2e", rel(r2$p, st[colnames(Ls), "p_model"])),
            sprintf("R %d / Python %d", nrow(pr), pp$n_pairs), sprintf("R %.0f / Python %.0f", sum(pr$diff), pp$sum_diff)))
print(res, right = FALSE)
data.table::fwrite(res, "tests/crosscheck_results.csv")

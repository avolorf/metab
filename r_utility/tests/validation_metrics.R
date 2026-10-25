#!/usr/bin/env Rscript
# Количественная оценка корректности на нескольких контрольных наборах (разные seed):
# чувствительность и доля ложных находок, восстановление факторов разведения и индекса загрязнения,
# объединение аддуктов, CV на случайных метках с отбором признаков внутри и вне CV.
#   Rscript tests/validation_metrics.R [число наборов]
if (!isTRUE(l10n_info()$`UTF-8`)) invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
options(urometabr.dir = "R")
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
source("run_analysis.R", encoding = "UTF-8"); source("tests/make_control_data.R", encoding = "UTF-8")
a <- commandArgs(trailingOnly = TRUE); nset <- if (length(a)) as.integer(a[1]) else 5
rows <- list()
for (k in seq_len(nset)) for (null in c(FALSE, TRUE)) {
  d <- file.path(tempdir(), sprintf("vm%d_%s", k, null)); unlink(d, recursive = TRUE)
  make_control_data(d, seed = 100 + k, null = null, series_file = "example/contaminant_series.csv")
  cfg <- read_config("config/config_example.yaml", quick = TRUE)
  cfg$input$data <- file.path(d, "feature_table.csv"); cfg$input$metadata <- file.path(d, "metadata.csv")
  cfg$output <- file.path(d, "out"); cfg$multivariate$run <- FALSE; cfg$report$html <- FALSE; cfg$report$fig_dpi <- 72
  invisible(capture.output(R <- run_pipeline(cfg)))
  tr <- read.csv(file.path(d, "truth.csv"), stringsAsFactors = FALSE); rownames(tr) <- tr$feature
  ts <- read.csv(file.path(d, "truth_samples.csv"), colClasses = c(sample = "character"))
  kept <- rownames(R$features)
  r <- R$res$BC_vs_Control; f <- rownames(r)[r$finding]; tt <- tr[f, "type"]
  truepos <- tr$feature[tr$type %in% c("BC_up", "BC_down")]
  # CV на случайных метках среди контролей
  ctrl <- R$meta$sample[R$meta$study & R$meta$group == "Control"]; set.seed(k); g1 <- sample(ctrl, 22); g0 <- setdiff(ctrl, g1)
  y <- c(rep(1, 22), rep(0, length(g0))); X <- R$Lc[c(g1, g0), ]
  cfg$multivariate$cv_repeats <- 5
  set.seed(k); auc_in <- mean(replicate(5, auc_score(y, cv_predict(X, y, "panel", 5, 5)$pred)))
  sel <- order(mw_p_fast(X, y))[1:5]
  set.seed(k); auc_out <- mean(replicate(5, auc_score(y, cv_predict(X[, sel], y, "covariate", 5, 5)$pred)))
  rows[[length(rows) + 1]] <- data.frame(
    set = k, null_effects = null,
    noisy_removed = mean(!tr$feature[tr$type == "noisy"] %in% kept),
    rare_blank_removed = mean(!tr$feature[tr$type %in% c("rare", "blank")] %in% kept),
    pqn_r = cor(log(R$meta[ts$sample, "dilution"]), log(ts$dilution_true)),
    contam_rho = cor(R$meta[ts$sample, "contam_index"], ts$contam_true, method = "spearman"),
    n_findings_BC = length(f), true_findings = sum(tt %in% c("BC_up", "BC_down")),
    sensitivity = if (null) NA else sum(f %in% truepos) / length(truepos),
    false_findings = sum(!tt %in% c("BC_up", "BC_down", "drug_BC")),
    age_only_q_MW_PC = sum(R$res$PC_vs_Control[intersect(rownames(R$res$PC_vs_Control), tr$feature[tr$type == "age_only"]), "q_MW"] < 0.05),
    age_only_q_model_PC = sum(R$res$PC_vs_Control[intersect(rownames(R$res$PC_vs_Control), tr$feature[tr$type == "age_only"]), "q_model"] < 0.05),
    findings_PC = sum(R$res$PC_vs_Control$finding),
    AUC_random_selection_inside_CV = auc_in, AUC_random_selection_outside_CV = auc_out)
  cat(".")
}
res <- do.call(rbind, rows)
data.table::fwrite(res, "tests/validation_metrics.csv")
print(res, digits = 3)

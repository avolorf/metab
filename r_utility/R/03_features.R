# Группировка ионов в соединения, индекс экзогенного загрязнения, категории экспозиций.

ADDUCTS <- c("[M+H]+" = 1.007276, "[M+Na]+" = 22.989218, "[M+K]+" = 38.963158, "[M+NH4]+" = 18.033823,
             "[M+H-H2O]+" = 1.007276 - 18.010565, "[M+H-2H2O]+" = 1.007276 - 2 * 18.010565)
RELATIONS <- c("13C" = 1.003355, "Na-H" = 21.981942, "K-H" = 37.955882, "NH4-H" = 17.026547, "H2O" = 18.010565,
               "2H2O" = 36.021130, "CO2" = 43.989829, "HCOOH" = 46.005479, "CO" = 27.994915)
ELEM <- c(C = 12, H = 1.00782503, N = 14.0030740, O = 15.9949146, P = 30.9737620, S = 31.9720707,
          Cl = 34.9688527, Br = 78.9183371, Na = 22.9897693, K = 38.9637064)
EXPOSURE <- c("contaminant", "drug", "tobacco", "food")

formula_mass <- function(f) {
  m <- regmatches(f, gregexpr("([A-Z][a-z]?)(\\d*)", f))[[1]]
  sum(vapply(m, function(t) {
    e <- sub("\\d+$", "", t); n <- sub("^[A-Za-z]+", "", t)
    ELEM[[e]] * ifelse(nzchar(n), as.numeric(n), 1)
  }, 0))
}
ppm_err <- function(obs, theo) (obs - theo) / theo * 1e6

# Соединения — связные компоненты графа «|ΔRT| ≤ rt_tol и ρ Спирмена ≥ min_corr».
# Представитель — признак с наибольшей медианной интенсивностью.
group_features <- function(L, info, raw_median, rt_tol = 0.02, min_corr = 0.8) {
  feats <- colnames(L); rt <- info[feats, "rt_min"]; mz <- info[feats, "mz"]
  n <- length(feats); o <- order(rt)
  R <- apply(L, 2, rank)
  R <- scale(R, center = TRUE, scale = apply(R, 2, function(v) sqrt(mean((v - mean(v))^2))) + 1e-12)
  parent <- seq_len(n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  edges <- list()
  for (a in seq_len(n)) {
    i <- o[a]; b <- a + 1
    while (b <= n && rt[o[b]] - rt[i] <= rt_tol) {
      j <- o[b]
      r <- sum(R[, i] * R[, j]) / nrow(R)
      if (r >= min_corr) {
        edges[[length(edges) + 1]] <- c(i, j)
        pi <- find(i); pj <- find(j)
        if (pi != pj) parent[pi] <- pj
      }
      b <- b + 1
    }
  }
  root <- vapply(seq_len(n), find, 1L)
  df <- data.frame(feature = feats, rt_min = rt, mz = mz, median_int = raw_median[feats], root = root,
                   stringsAsFactors = FALSE)
  ordint <- order(-df$median_int)
  rep <- df$feature[ordint][!duplicated(df$root[ordint])]
  names(rep) <- df$root[ordint][!duplicated(df$root[ordint])]
  df$representative <- rep[as.character(df$root)]
  roots <- unique(df$root[order(df$rt_min)])
  df$compound <- sprintf("C%05d", match(df$root, roots))
  df$n_ions <- as.integer(table(df$compound)[df$compound])
  rel <- vector("list", n)
  for (e in edges) {
    d <- abs(mz[e[1]] - mz[e[2]])
    for (nm in names(RELATIONS)) if (abs(d - RELATIONS[[nm]]) < 0.003) {
      lo <- e[which.min(mz[e])]; hi <- e[which.max(mz[e])]
      rel[[hi]] <- c(rel[[hi]], sprintf("%s от %.4f", nm, mz[lo]))
      rel[[lo]] <- c(rel[[lo]], sprintf("%.4f (+%s)", mz[hi], nm))
    }
  }
  df$ion_relations <- vapply(rel, function(v) paste(v, collapse = ";"), "")
  df$is_representative <- df$feature == df$representative
  rownames(df) <- df$feature
  df$root <- NULL
  df
}

# Члены гомологического ряда экзогенных соединений (файл: label, formula, ion, rt_min, in_index)
match_series <- function(info, path, ppm = 5) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  s <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8")
  rows <- list()
  for (k in seq_len(nrow(s))) {
    theo <- formula_mass(s$formula[k]) + ADDUCTS[[s$ion[k]]]
    hit <- which(abs(ppm_err(info$mz, theo)) <= ppm & abs(info$rt_min - s$rt_min[k]) <= 0.05)
    for (h in hit) rows[[length(rows) + 1]] <- data.frame(
      feature = info$feature[h], label = s$label[k], formula = s$formula[k], ion = s$ion[k],
      rt_min = info$rt_min[h], mz = info$mz[h], ppm = round(ppm_err(info$mz[h], theo), 2),
      in_index = as.integer(s$in_index[k]), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out <- out[!duplicated(out$feature), ]
  rownames(out) <- out$feature
  out
}

# Индекс загрязнения: средний z-score членов ряда (z по исследуемым пробам); high/low по медиане
contamination_index <- function(L, series, study_ids) {
  f <- intersect(series$feature[series$in_index == 1], colnames(L))
  if (!length(f)) return(NULL)
  Z <- L[, f, drop = FALSE]
  mu <- colMeans(Z[study_ids, , drop = FALSE]); sdv <- apply(Z[study_ids, , drop = FALSE], 2, stats::sd)
  Z <- sweep(sweep(Z, 2, mu, "-"), 2, sdv, "/")
  idx <- rowMeans(Z)
  thr <- stats::median(idx[study_ids])
  list(index = idx, level = ifelse(idx > thr, "high", "low"), threshold = thr, features = f)
}

exposure_category <- function(info, path) {
  cat <- stats::setNames(rep("", nrow(info)), info$feature)
  if (is.null(path) || !file.exists(path)) return(cat)
  kw <- data.table::fread(path, data.table = FALSE, encoding = "UTF-8")
  for (k in seq_len(nrow(kw))) {
    m <- grepl(kw$pattern[k], info$annotation, ignore.case = TRUE, perl = TRUE) & cat == ""
    cat[m] <- kw$category[k]
  }
  cat
}

# Категория признака и соединения: contaminant > drug > tobacco > food > прочие
categorize <- function(info, groups, series, kw_cat) {
  cat <- stats::setNames(rep("unknown", nrow(info)), info$feature)
  cat[info$annotation_conf == "high"] <- "annotated (high score)"
  cat[info$annotation_conf == "low"] <- "putative (low score)"
  hit <- kw_cat != ""
  cat[hit] <- ifelse(info$annotation_conf[hit] == "high", kw_cat[hit], paste0("possible ", kw_cat[hit], " (low score)"))
  if (!is.null(series)) cat[intersect(series$feature, names(cat))] <- "contaminant"
  prio <- c("contaminant", "drug", "tobacco", "food", "possible drug (low score)", "possible tobacco (low score)",
            "possible food (low score)", "endogenous", "annotated (high score)", "putative (low score)", "unknown")
  g <- groups[, c("feature", "compound")]
  g$category <- cat[g$feature]
  g$rank <- match(g$category, prio)
  best <- g[order(g$rank), ]; best <- best[!duplicated(best$compound), ]
  g$compound_category <- stats::setNames(best$category, best$compound)[g$compound]
  rownames(g) <- g$feature
  g[, c("compound", "category", "compound_category")]
}

panel_excluded <- function(cat) cat %in% EXPOSURE | startsWith(cat, "possible drug") | startsWith(cat, "possible tobacco")

display_name <- function(info) {
  tag <- sprintf("%.4f@%.2f", info$mz, info$rt_min)
  stats::setNames(ifelse(info$annotation_conf == "high", paste(substr(info$annotation, 1, 32), tag),
                         paste0("unknown ", tag, ifelse(info$annotation_conf == "low", "*", ""))), info$feature)
}

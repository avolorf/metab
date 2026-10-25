# Модуль 1. Оптимизация параметров (IPO) и таблица пиков (xcms: centWave → Obiwarp → PeakDensity →
# fillChromPeaks) по сырым данным. Исправлено относительно этапа 1: пути из настроек, один режим
# параллельных вычислений, minFraction — параметр (по умолчанию 0,5), группы образцов из метаданных.
# Нужны Bioconductor-пакеты xcms и MSnbase (IPO — только при xcms$ipo: true).

bp_param <- function(cores) {
  if (cores <= 1) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") BiocParallel::SnowParam(cores) else BiocParallel::MulticoreParam(cores)
}

module1_peaktable <- function(s) {
  for (p in c("xcms", "MSnbase", "BiocParallel"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Для обработки сырых данных нужен пакет ", p, ": BiocManager::install(c('xcms', 'MSnbase', 'IPO'))")
  out <- module_dir(s, 1); x <- s$xcms
  inp <- load_module(s, 0); meta <- inp$meta; files <- inp$files[meta$file]
  bp <- bp_param(s$project$cores)
  BiocParallel::register(bp)
  P <- list(ppm = x$ppm, min_peakwidth = x$min_peakwidth, max_peakwidth = x$max_peakwidth, snthresh = x$snthresh,
            noise = x$noise, prefilter = x$prefilter_k, value_of_prefilter = x$prefilter_i, mzdiff = -0.001,
            integrate = 1L, fitgauss = FALSE, mzCenterFun = "wMean")
  R <- list(binsize = x$obiwarp_binsize, bw = x$bw, mzwid = x$mzwid)
  if (isTRUE(x$ipo)) {
    if (!requireNamespace("IPO", quietly = TRUE)) stop("xcms$ipo: true, но пакет IPO не установлен")
    qc <- files[meta$kind == "QC"]
    if (length(qc) < 2) stop("Для оптимизации IPO нужны минимум 2 файла QC")
    log_msg("IPO: оптимизация centWave по ", length(qc), " файлам QC")
    pp <- IPO::getDefaultXcmsSetStartingParams("centWave")
    pp$ppm <- c(max(1, x$ppm - 2), x$ppm + 5); pp$min_peakwidth <- c(3, 15); pp$max_peakwidth <- c(20, 40)
    pp$noise <- c(100, 1000); pp$value_of_prefilter <- c(3, 800)
    rp <- IPO::optimizeXcmsSet(files = qc, params = pp, nSlaves = s$project$cores, subdir = NULL, plot = FALSE)
    rg <- IPO::optimizeRetGroup(xset = rp$best_settings$xset, params = IPO::getDefaultRetGroupStartingParams(),
                                nSlaves = 1, subdir = NULL, plot = FALSE)
    P <- utils::modifyList(P, rp$best_settings$parameters)
    R <- list(binsize = rg$best_settings$profStep, bw = rg$best_settings$bw, mzwid = rg$best_settings$mzwid)
    saveRDS(list(peakpicking = rp$best_settings$parameters, retgroup = rg$best_settings), file.path(out, "ipo_params.rds"))
  }
  pd <- data.frame(sample_name = meta$file, group = meta$group, stringsAsFactors = FALSE)
  raw <- MSnbase::readMSData(files, pdata = methods::new("NAnnotatedDataFrame", pd), mode = "onDisk")
  cwp <- xcms::CentWaveParam(ppm = P$ppm, peakwidth = c(P$min_peakwidth, P$max_peakwidth), snthresh = P$snthresh,
                             prefilter = c(P$prefilter, P$value_of_prefilter), noise = P$noise, mzdiff = P$mzdiff,
                             integrate = as.integer(P$integrate), fitgauss = as.logical(P$fitgauss), mzCenterFun = P$mzCenterFun)
  log_msg("Детектирование пиков (centWave)…")
  xd <- xcms::findChromPeaks(raw, param = cwp, BPPARAM = bp)
  log_msg("Выравнивание времён удерживания (Obiwarp)…")
  xd <- xcms::adjustRtime(xd, param = xcms::ObiwarpParam(binSize = R$binsize))
  xd <- xcms::groupChromPeaks(xd, param = xcms::PeakDensityParam(sampleGroups = pd$group, bw = R$bw,
                                                                  minFraction = x$min_fraction, binSize = R$mzwid))
  xd <- xcms::fillChromPeaks(xd, param = xcms::ChromPeakAreaParam(), BPPARAM = bp)
  ft <- xcms::featureValues(xd, value = "into"); fd <- xcms::featureDefinitions(xd)
  X <- t(ft); rownames(X) <- meta$file; X[!is.finite(X) | X <= 0] <- NA
  info <- data.frame(id = sprintf("F%05d", seq_len(ncol(X))), original = rownames(fd), mz = fd$mzmed,
                     rt_min = fd$rtmed / 60, name = "", stringsAsFactors = FALSE)
  colnames(X) <- info$id
  log_msg(sprintf("Таблица пиков: %d признаков, пропусков %.1f%% (minFraction = %g)", ncol(X), 100 * mean(is.na(X)), x$min_fraction))
  write_csv(cbind(file = rownames(X), as.data.frame(X)), file.path(out, "peak_table.csv"))
  write_csv(info, file.path(out, "features.csv"))
  saveRDS(xd, file.path(out, "xcms_object.rds"))
  save_module(s, 1, list(X = X, info = info, meta = meta, type = "raw"))
}

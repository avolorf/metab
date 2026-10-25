# Необязательный модуль: таблица признаков из сырых данных (mzML/mzXML/CDF) средствами xcms.
# Исправленная версия скрипта макета 1 этапа (см. пояснительную записку, раздел 2.3).
# Требует Bioconductor: BiocManager::install(c("xcms", "MSnbase", "IPO"));
# в основном сценарии утилиты не используется (таблица экспортируется из MS-DIAL).

xcms_peak_table <- function(cfg) {
  for (p in c("xcms", "MSnbase", "BiocParallel"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Для input$format = 'xcms' установите пакет ", p)
  d <- cfg$input$mzml_dir
  files <- list.files(d, pattern = "\\.(mzML|mzXML|CDF)$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("В ", d, " нет файлов mzML/mzXML/CDF")
  x <- cfg$xcms %||% list()
  # одна схема параллельных вычислений (в макете SnowParam сразу перезаписывался SerialParam)
  bp <- if ((cfg$cores %||% 1) > 1) BiocParallel::SnowParam(cfg$cores) else BiocParallel::SerialParam()
  pd <- data.frame(sample_name = tools::file_path_sans_ext(basename(files)), stringsAsFactors = FALSE)
  raw <- MSnbase::readMSData(files, pdata = methods::new("NAnnotatedDataFrame", pd), mode = "onDisk")
  cwp <- xcms::CentWaveParam(ppm = x$ppm %||% 5, peakwidth = unlist(x$peakwidth %||% c(5, 30)),
                             snthresh = x$snthresh %||% 10, prefilter = unlist(x$prefilter %||% c(3, 1000)),
                             noise = x$noise %||% 1000, mzdiff = x$mzdiff %||% 0.001, integrate = 1L)
  xd <- xcms::findChromPeaks(raw, param = cwp, BPPARAM = bp)
  xd <- xcms::adjustRtime(xd, param = xcms::ObiwarpParam(binSize = x$obiwarp_binsize %||% 0.5))
  # minFraction задаётся явно (0.5), а не берётся из IPO; sampleGroups — по типу пробы
  grp <- ifelse(grepl(cfg$design$qc_pattern, pd$sample_name, ignore.case = TRUE), "QC", "sample")
  pdp <- xcms::PeakDensityParam(sampleGroups = grp, bw = x$bw %||% 5, minFraction = x$min_fraction %||% 0.5,
                                binSize = x$mzwid %||% 0.015)
  xd <- xcms::groupChromPeaks(xd, param = pdp)
  xd <- xcms::fillChromPeaks(xd, param = xcms::ChromPeakAreaParam(), BPPARAM = bp)
  ft <- xcms::featureValues(xd, value = "into")
  fd <- xcms::featureDefinitions(xd)
  X <- t(ft); X[X == 0] <- NA
  rownames(X) <- pd$sample_name
  fid <- make.unique(sprintf("%.3f_%.5f_Unknown", fd$rtmed / 60, fd$mzmed))
  colnames(X) <- fid
  saveRDS(xd, file.path(cfg$output, "xcms_object.rds"))
  list(X = X, info = feature_info(fid, fd$rtmed / 60, fd$mzmed, rep("Unknown", length(fid))),
       kind = classify_injections(rownames(X), rep("", nrow(X)), cfg), inj_order = NULL, n_dup = 0L)
}

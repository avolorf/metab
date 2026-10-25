# Формирование отчёта: таблицы (CSV, XLSX), HTML-отчёт с рисунками, журнал, сведения о сессии,
# контрольные суммы результатов.

write_tables <- function(sheets, out) {
  dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)
  wb <- openxlsx::createWorkbook()
  for (k in names(sheets)) {
    d <- sheets[[k]]
    if (is.null(d) || !nrow(d)) next
    if (!is.null(rownames(d)) && !identical(rownames(d), as.character(seq_len(nrow(d)))) && !"id" %in% names(d))
      d <- cbind(id = rownames(d), d)
    data.table::fwrite(d, file.path(out, "tables", paste0(k, ".csv")))
    sh <- substr(k, 1, 31)
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, d)
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
  }
  openxlsx::saveWorkbook(wb, file.path(out, "results.xlsx"), overwrite = TRUE)
}

html_table <- function(d, digits = 3) {
  if (is.null(d) || !nrow(d)) return("<p><i>нет данных</i></p>")
  f <- function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "g") else htmlesc(as.character(v))
  body <- apply(as.data.frame(lapply(d, f), stringsAsFactors = FALSE), 1,
                function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", paste0("<th>", htmlesc(names(d)), "</th>", collapse = ""), "</tr></thead><tbody>",
         paste(body, collapse = "\n"), "</tbody></table>")
}

htmlesc <- function(s) { s <- gsub("&", "&amp;", s); s <- gsub("<", "&lt;", s); gsub(">", "&gt;", s) }

img_tag <- function(png, caption) {
  if (is.null(png) || !file.exists(png)) return("")
  b64 <- jsonlite::base64_enc(readBin(png, "raw", file.info(png)$size))
  sprintf("<figure><img src=\"data:image/png;base64,%s\"/><figcaption>%s</figcaption></figure>", b64, htmlesc(caption))
}

write_html_report <- function(out, cfg, sections) {
  css <- "body{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;color:#222}
table{border-collapse:collapse;font-size:12px;margin:8px 0}td,th{border:1px solid #ccc;padding:3px 6px}
th{background:#f0f0f0}img{max-width:100%}figcaption{font-size:12px;color:#555}pre{background:#f7f7f7;padding:8px;font-size:11px;white-space:pre-wrap}"
  body <- vapply(sections, function(s) paste0("<h2>", htmlesc(s$title), "</h2>", s$html), "")
  html <- paste0("<!doctype html><html lang=\"ru\"><head><meta charset=\"utf-8\"><title>Отчёт UroMetabR</title><style>", css,
                 "</style></head><body><h1>Отчёт об обработке данных ВЭЖХ-МС (UroMetabR ", UROMETABR_VERSION, ")</h1>",
                 "<p>Дата: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "; данные: ", htmlesc(basename(cfg$input$data %||% "")),
                 "; seed: ", cfg$seed, "</p>", paste(body, collapse = "\n"), "</body></html>")
  writeLines(html, file.path(out, "report.html"), useBytes = TRUE)
}

write_manifest <- function(out) {
  f <- sort(list.files(file.path(out, "tables"), full.names = TRUE))
  m <- data.frame(file = basename(f), md5 = unname(md5_files(f)))
  data.table::fwrite(m, file.path(out, "manifest_md5.csv"))
  m
}

write_session_info <- function(out) {
  si <- utils::capture.output(utils::sessionInfo())
  writeLines(c(paste("UroMetabR", UROMETABR_VERSION), si), file.path(out, "session_info.txt"))
}

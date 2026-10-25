# Графический интерфейс UroMetabR (Shiny). Запуск: Start_UroMetabR.bat (Windows), start_urometabr.sh
# (Linux/macOS) или в R: shiny::runApp("app"). Пользователь выбирает папку проекта, файлы данных и
# параметры; обработка идёт в отдельном процессе R, журнал обновляется на экране.

library(shiny)
root <- getOption("urometabr.root", normalizePath(file.path(getwd(), if (basename(getwd()) == "app") ".." else ".")))
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8"); load_urometabr(root)

is_win <- .Platform$OS.type == "windows"
pick_dir <- function() if (is_win) utils::choose.dir() else NA
pick_file <- function() if (is_win) tryCatch(utils::choose.files(multi = FALSE), error = function(e) NA) else NA
num <- function(id, label, value, step = NULL, min = NA, max = NA) numericInput(id, label, value, min = min, max = max, step = step, width = "100%")

ui <- fluidPage(
  tags$head(tags$script(HTML("$(document).on('shiny:value', function(e){ if (e.name === 'log_text') setTimeout(function(){ var d = document.getElementById('log_box'); d.scrollTop = d.scrollHeight; }, 50); });")),
            tags$style("body{font-size:14px} #fig_img img{max-width:900px} .mono{font-family:Consolas,monospace;font-size:12px;white-space:pre-wrap;background:#f6f6f6;padding:8px;max-height:420px;overflow-y:auto}")),
  titlePanel(div(h3(PROGRAM_NAME), tags$small(paste("UroMetabR", UROMETABR_VERSION, "от", UROMETABR_DATE)))),
  fluidRow(column(7, textInput("proj", "Папка проекта (результаты всех модулей сохраняются в ней)", width = "100%")),
           column(2, br(), actionButton("proj_pick", "Выбрать…", width = "100%")),
           column(3, br(), actionButton("proj_open", "Открыть / создать проект", class = "btn-primary", width = "100%"))),
  textOutput("proj_status"),
  tabsetPanel(id = "tabs",
    tabPanel("1. Данные", br(),
      radioButtons("in_type", "Тип входных данных", c("Готовая таблица признаков (CSV, TSV, TXT, экспорт MS-DIAL)" = "table",
                                                      "Сырые данные (папка с mzML / mzXML / CDF; нужен пакет xcms)" = "raw")),
      conditionalPanel("input.in_type == 'table'",
        fluidRow(column(8, textInput("table_file", "Файл таблицы признаков", width = "100%")),
                 column(2, br(), actionButton("table_pick", "Выбрать…", width = "100%")),
                 column(2, fileInput("table_up", "или загрузить", buttonLabel = "Файл…", placeholder = "")))),
      conditionalPanel("input.in_type == 'raw'",
        fluidRow(column(8, textInput("raw_dir", "Папка с файлами mzML / mzXML / CDF", width = "100%")),
                 column(2, br(), actionButton("raw_pick", "Выбрать…", width = "100%")))),
      fluidRow(column(8, textInput("meta_file", "Файл метаданных (столбцы: file, group; необязательно: sample, injection_order, replicate, include, sex, age)", width = "100%")),
               column(2, br(), actionButton("meta_pick", "Выбрать…", width = "100%")),
               column(2, fileInput("meta_up", "или загрузить", buttonLabel = "Файл…", placeholder = ""))),
      actionButton("check", "Проверить данные (модуль 0)", class = "btn-info"), br(), br(),
      tableOutput("check_tab")),
    tabPanel("2. Параметры", br(),
      fluidRow(
        column(4, h4("Фильтр (модуль 3)"),
               num("f_rsd", "Макс. RSD в QC", 0.30, 0.05, 0, 1), num("f_grp", "Мин. доля присутствия в группе", 0.5, 0.05, 0, 1),
               num("f_qc", "Мин. доля обнаружения в QC", 0.8, 0.05, 0, 1), num("f_blank", "Мин. отношение к бланку", 5, 1, 1),
               h4("Пропуски и коррекция (5–6)"),
               selectInput("imp", "Заполнение пропусков", c("Случайный лес, алгоритм missForest" = "rf", "Пакет missForest" = "missForest", "½ минимума" = "halfmin")),
               num("imp_iter", "Макс. итераций", 10, 1, 1), num("imp_tree", "Деревьев", 100, 10, 10),
               selectInput("norm", "Нормализация", c("нет" = "none", "PQN" = "PQN")),
               selectInput("corr", "Коррекция сигнала", c("EigenMS" = "eigenms", "нет" = "none")), num("eig_perm", "Перестановок EigenMS", 200, 50, 20)),
        column(4, h4("Отбор признаков (модуль 7)"),
               selectInput("padj", "Поправка на множественные сравнения", c("Холм–Бонферрони" = "holm", "Бенджамини–Хохберг (FDR)" = "BH",
                                                                            "без поправки" = "none", "как на этапе 1 (без поправки по признакам)" = "stage1")),
               num("alpha", "Уровень значимости", 0.05, 0.01, 0.001, 0.2),
               checkboxGroupInput("models", "Модели", c("KNN" = "knn", "SVM (радиальное ядро)" = "svmRadial", "PLS" = "pls", "Случайный лес" = "rf"),
                                  c("knn", "svmRadial", "pls", "rf"), inline = TRUE),
               num("top_n", "Важнейших признаков на модель", 50, 5, 5), num("min_freq", "Мин. число моделей", 2, 1, 1, 4),
               selectInput("rfe", "РОП: классификатор", c("наивный Байес" = "nb", "случайный лес" = "rf")),
               checkboxInput("nested", "Вложенная кросс-валидация (несмещённая оценка точности)", TRUE),
               num("n_folds", "Блоков вложенной CV", 5, 1, 3, 10), num("n_rep", "Повторов вложенной CV", 2, 1, 1, 10)),
        column(4, h4("Общие"),
               textInput("pname", "Название проекта", "Проект"), num("seed", "Seed (одно число для всех модулей)", 20261026, 1),
               num("cores", "Потоков", 1, 1, 1, 64), num("dpi", "Разрешение рисунков, dpi", 600, 50, 72, 1200),
               checkboxGroupInput("fmt", "Форматы рисунков", c("PNG" = "png", "JPEG" = "jpeg"), c("png", "jpeg"), inline = TRUE),
               textInput("db", "Папка базы признаков (дозапись)", "features_db"),
               h4("xcms (модуль 1, сырые данные)"),
               num("x_ppm", "ppm", 5, 1, 1), num("x_pw1", "Мин. ширина пика, с", 5, 1, 1), num("x_pw2", "Макс. ширина пика, с", 30, 1, 2),
               num("x_minfrac", "minFraction", 0.5, 0.05, 0, 1), checkboxInput("x_ipo", "Оптимизировать параметры (IPO)", FALSE))),
      actionButton("save", "Сохранить настройки", class = "btn-primary"), textOutput("save_status"), br()),
    tabPanel("3. Запуск", br(),
      fluidRow(column(3, selectInput("from", "С модуля", setNames(0:9, paste(0:9, MODULES$title)), 0)),
               column(3, selectInput("to", "по модуль", setNames(0:9, paste(0:9, MODULES$title)), 9)),
               column(3, br(), actionButton("run", "Запустить обработку", class = "btn-success", width = "100%"))),
      h4(textOutput("run_state")), tableOutput("mod_tab"), div(id = "log_box", class = "mono", textOutput("log_text"))),
    tabPanel("4. Результаты", br(),
      actionButton("open_report", "Открыть отчёт (HTML)", class = "btn-primary"), actionButton("open_dir", "Открыть папку проекта"), br(), br(),
      h4("Число признаков и средняя AUROC по этапам"), tableOutput("res_steps"),
      h4("Вложенная кросс-валидация"), tableOutput("res_nested"),
      h4("Отобранные признаки (РОП)"), tableOutput("res_rfe"),
      selectInput("fig", "Рисунок", character(0), width = "60%"), imageOutput("fig_img", height = "auto")),
    tabPanel("Справка", br(), includeMarkdown(file.path(root, "docs", "user_guide.md"))))
)

server <- function(input, output, session) {
  rv <- reactiveValues(s = NULL, running = FALSE, started = NULL, msg = "")
  observeEvent(input$proj_pick, { d <- pick_dir(); if (!is.na(d)) updateTextInput(session, "proj", value = d)
    else showNotification("Выбор папки окном доступен в Windows; введите путь вручную", type = "warning") })
  for (k in c("table", "meta")) local({ kk <- k
    observeEvent(input[[paste0(kk, "_pick")]], { f <- pick_file(); if (length(f) && !is.na(f)) updateTextInput(session, paste0(kk, "_file"), value = f)
      else showNotification("Выбор файла окном доступен в Windows; введите путь или загрузите файл", type = "warning") })
    observeEvent(input[[paste0(kk, "_up")]], {
      req(rv$s); up <- input[[paste0(kk, "_up")]]; d <- file.path(rv$s$project$dir, "input"); dir.create(d, showWarnings = FALSE)
      f <- file.path(d, up$name); file.copy(up$datapath, f, overwrite = TRUE); updateTextInput(session, paste0(kk, "_file"), value = f) }) })
  observeEvent(input$raw_pick, { d <- pick_dir(); if (!is.na(d)) updateTextInput(session, "raw_dir", value = d) })

  fill_ui <- function(s) {
    upd <- function(id, v) updateNumericInput(session, id, value = v)
    updateTextInput(session, "pname", value = s$project$name); upd("seed", s$project$seed); upd("cores", s$project$cores)
    updateRadioButtons(session, "in_type", selected = s$input$type)
    updateTextInput(session, "table_file", value = s$input$table_file); updateTextInput(session, "meta_file", value = s$input$metadata_file)
    updateTextInput(session, "raw_dir", value = s$input$raw_dir)
    upd("f_rsd", s$filter$max_qc_rsd); upd("f_grp", s$filter$min_group_detect); upd("f_qc", s$filter$min_qc_detect); upd("f_blank", s$filter$blank_ratio)
    updateSelectInput(session, "imp", selected = s$impute$method); upd("imp_iter", s$impute$maxiter); upd("imp_tree", s$impute$ntree)
    updateSelectInput(session, "norm", selected = s$correct$normalization); updateSelectInput(session, "corr", selected = s$correct$method); upd("eig_perm", s$correct$eigen_perm)
    updateSelectInput(session, "padj", selected = s$select$p_adjust); upd("alpha", s$select$alpha)
    updateCheckboxGroupInput(session, "models", selected = s$select$models); upd("top_n", s$select$top_n); upd("min_freq", s$select$min_freq)
    updateSelectInput(session, "rfe", selected = s$select$rfe_method); updateCheckboxInput(session, "nested", value = isTRUE(s$select$nested_cv))
    upd("n_folds", s$select$nested_folds); upd("n_rep", s$select$nested_repeats); upd("dpi", s$visual$dpi)
    updateCheckboxGroupInput(session, "fmt", selected = s$visual$formats); updateTextInput(session, "db", value = s$report$database)
    upd("x_ppm", s$xcms$ppm); upd("x_pw1", s$xcms$min_peakwidth); upd("x_pw2", s$xcms$max_peakwidth); upd("x_minfrac", s$xcms$min_fraction)
    updateCheckboxInput(session, "x_ipo", value = isTRUE(s$xcms$ipo))
  }
  collect <- function() {
    s <- rv$s
    s$project$name <- input$pname; s$project$seed <- input$seed; s$project$cores <- input$cores
    s$input$type <- input$in_type; s$input$table_file <- input$table_file; s$input$metadata_file <- input$meta_file; s$input$raw_dir <- input$raw_dir
    s$filter$max_qc_rsd <- input$f_rsd; s$filter$min_group_detect <- input$f_grp; s$filter$min_qc_detect <- input$f_qc; s$filter$blank_ratio <- input$f_blank
    s$impute$method <- input$imp; s$impute$maxiter <- input$imp_iter; s$impute$ntree <- input$imp_tree
    s$correct$normalization <- input$norm; s$correct$method <- input$corr; s$correct$eigen_perm <- input$eig_perm
    s$select$p_adjust <- input$padj; s$select$alpha <- input$alpha; s$select$models <- input$models; s$select$top_n <- input$top_n
    s$select$min_freq <- input$min_freq; s$select$rfe_method <- input$rfe; s$select$nested_cv <- input$nested
    s$select$nested_folds <- input$n_folds; s$select$nested_repeats <- input$n_rep
    s$visual$dpi <- input$dpi; s$visual$formats <- input$fmt; s$report$database <- input$db
    s$xcms$ppm <- input$x_ppm; s$xcms$min_peakwidth <- input$x_pw1; s$xcms$max_peakwidth <- input$x_pw2
    s$xcms$min_fraction <- input$x_minfrac; s$xcms$ipo <- input$x_ipo
    s$run$from <- as.integer(input$from); s$run$to <- as.integer(input$to)
    s
  }
  save_s <- function() { s <- collect(); check_settings(s); write_settings(s, s$settings_file); rv$s <- read_settings(s$settings_file); rv$s }

  observeEvent(input$proj_open, {
    d <- trimws(input$proj); if (!nzchar(d)) { showNotification("Укажите папку проекта", type = "error"); return() }
    dir.create(d, showWarnings = FALSE, recursive = TRUE); f <- file.path(d, "settings.yaml")
    if (!file.exists(f)) write_settings(default_settings(), f)
    rv$s <- tryCatch(read_settings(f), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    if (!is.null(rv$s)) fill_ui(rv$s)
  })
  output$proj_status <- renderText(if (is.null(rv$s)) "Проект не открыт" else paste("Открыт проект:", rv$s$project$dir, "— настройки:", rv$s$settings_file))
  observeEvent(input$save, { req(rv$s); r <- tryCatch({ save_s(); "Настройки сохранены" }, error = function(e) conditionMessage(e)); output$save_status <- renderText(r) })
  observeEvent(input$check, {
    req(rv$s); s <- tryCatch(save_s(), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }); req(s)
    log_start(s); r <- tryCatch({ module0_input(s); NULL }, error = function(e) conditionMessage(e))
    if (!is.null(r)) { showNotification(r, type = "error", duration = NULL); output$check_tab <- renderTable(data.frame(Ошибка = r)) }
    else output$check_tab <- renderTable(utils::read.csv(file.path(s$project$dir, "00_input", "input_check.csv"), encoding = "UTF-8", check.names = FALSE))
  })
  observeEvent(input$run, {
    req(rv$s); if (rv$running) return()
    s <- tryCatch(save_s(), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }); req(s)
    rs <- file.path(R.home("bin"), if (is_win) "Rscript.exe" else "Rscript")
    rv$started <- Sys.time(); rv$running <- TRUE
    system2(rs, c(shQuote(file.path(root, "run_analysis.R")), shQuote(s$settings_file), paste0("--from=", input$from), paste0("--to=", input$to)),
            wait = FALSE, stdout = file.path(s$project$dir, "logs", "console.txt"), stderr = file.path(s$project$dir, "logs", "console.txt"))
    updateTabsetPanel(session, "tabs", "3. Запуск")
  })
  latest <- function(pat) { if (is.null(rv$s)) return(NULL); f <- list.files(file.path(rv$s$project$dir, "logs"), pattern = pat, full.names = TRUE)
    if (!length(f)) return(NULL); f <- f[order(file.info(f)$mtime)]; if (!is.null(rv$started)) f <- f[file.info(f)$mtime >= rv$started - 1]
    if (length(f)) utils::tail(f, 1) else NULL }
  log_now <- reactivePoll(1000, session, function() { f <- latest("^run_.*\\.log$"); if (is.null(f)) 0 else file.info(f)$mtime },
                          function() { f <- latest("^run_.*\\.log$"); if (is.null(f)) character(0) else readLines(f, encoding = "UTF-8", warn = FALSE) })
  observe({ l <- log_now(); if (rv$running && any(grepl("UROMETABR_DONE|Обработка остановлена", l))) rv$running <- FALSE })
  output$run_state <- renderText({ l <- log_now()
    if (rv$running) "Идёт обработка…" else if (any(grepl("Обработка остановлена", l))) "Остановлено из-за ошибки — см. журнал"
    else if (any(grepl("UROMETABR_DONE", l))) "Обработка завершена" else "" })
  output$log_text <- renderText(paste(utils::tail(setdiff(log_now(), "UROMETABR_DONE"), 200), collapse = "\n"))
  output$mod_tab <- renderTable({ log_now(); f <- latest("^modules_.*\\.csv$"); if (is.null(f)) return(NULL)
    d <- utils::read.csv(f, encoding = "UTF-8", check.names = FALSE)
    stats::setNames(d, c("Модуль", "Название", "Статус", "Время, с", "Пиковая память R, МБ")) }, na = "")
  rd <- function(id, f) { log_now(); req(rv$s); p <- file.path(rv$s$project$dir, MODULES$dir[MODULES$id == id], f)
    if (file.exists(p)) utils::read.csv(p, encoding = "UTF-8", check.names = FALSE) }
  output$res_steps <- renderTable(rd(7, "auroc_steps.csv"), digits = 3)
  output$res_nested <- renderTable(rd(7, "nested_cv.csv"), digits = 3)
  output$res_rfe <- renderTable({ u <- rd(7, "uvf_stats.csv"); r <- rd(7, "ds_rfe.csv"); if (!is.null(u) && !is.null(r)) u[u$id %in% names(r), c("id", "mz", "rt_min", "name", "test", "p", "p_adj")] }, digits = 4, na = "")
  observe({ log_now(); req(rv$s); f <- list.files(file.path(rv$s$project$dir, c("02_qc", "06_correct", "08_visual")), pattern = "\\.png$", full.names = TRUE)
    updateSelectInput(session, "fig", choices = stats::setNames(f, basename(f))) })
  output$fig_img <- renderImage({ req(input$fig); list(src = input$fig, width = "100%") }, deleteFile = FALSE)
  observeEvent(input$open_report, { req(rv$s); f <- file.path(rv$s$project$dir, "09_report", "report.html")
    if (file.exists(f)) utils::browseURL(f) else showNotification("Отчёт ещё не сформирован", type = "warning") })
  observeEvent(input$open_dir, { req(rv$s); if (is_win) shell.exec(rv$s$project$dir) else utils::browseURL(rv$s$project$dir) })
}

shinyApp(ui, server)

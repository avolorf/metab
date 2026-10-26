# Запуск графического интерфейса: проверка/установка пакетов и открытие окна в браузере.
if (!isTRUE(l10n_info()$`UTF-8`) && .Platform$OS.type != "windows")      # кириллица в модулях: UTF-8
  for (loc in c("C.UTF-8", "en_US.UTF-8", "ru_RU.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(normalizePath(self)), ".."))
source(file.path(root, "R", "zz_setup.R"), encoding = "UTF-8")
if (length(missing_packages())) install_missing()
load_urometabr(root)
port <- as.integer(Sys.getenv("UROMETABR_PORT", "0")); if (!port) port <- NULL
message("\nUroMetabR запускается. Окно программы откроется в браузере.\n",
        "НЕ ЗАКРЫВАЙТЕ это чёрное окно, пока работаете с программой: его закрытие завершает работу UroMetabR.\n",
        "Если браузер не открылся, скопируйте адрес http://127.0.0.1:… из строки «Listening on» ниже в браузер.\n")
shiny::runApp(file.path(root, "app"), launch.browser = Sys.getenv("UROMETABR_NO_BROWSER") == "", port = port)

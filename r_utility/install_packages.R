# Установка пакетов, необходимых UroMetabR (однократно; нужен доступ к CRAN).
pk <- c("data.table", "yaml", "ggplot2", "patchwork", "ggrepel", "pheatmap", "openxlsx", "pls", "clue",
        "jsonlite", "digest", "testthat", "pROC")
miss <- pk[!vapply(pk, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
# Необязательно — обработка сырых данных mzML/mzXML/CDF (input$format: xcms):
# install.packages("BiocManager"); BiocManager::install(c("xcms", "MSnbase"))
cat("Готово. Не установлены:", paste(pk[!vapply(pk, requireNamespace, TRUE, quietly = TRUE)], collapse = ", "), "\n")

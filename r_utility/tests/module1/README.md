# Проверка модуля 1 (xcms) на открытых наборах faahKO и msdata

**Одна команда** (из папки утилиты `r_utility`, нужен интернет для установки пакетов Bioconductor):

    Rscript tests/module1_check.R --out=module1_check

В Docker-образе Bioconductor (из папки, где лежит `r_utility`):

    docker run --rm -v "$PWD":/work -w /work/r_utility bioconductor/bioconductor_docker:RELEASE_3_22 Rscript tests/module1_check.R --out=/work/module1_check

Скрипт ставит xcms, MSnbase, faahKO, msdata и недостающие пакеты CRAN; копирует 12 файлов faahKO (netCDF, 6 KO и 6 WT),
конвертирует их в mzML и mzXML, прогоняет утилиту на трёх форматах (CDF и mzXML — модули 0–3, mzML — модули 0–9) и проверяет
чтение файлов mzML прибора Sciex из msdata. Метаданные: `metadata_faahKO.csv` (file, group, injection_order),
`metadata_msdata_sciex.csv` (msdata используется только для проверки чтения формата: двух файлов недостаточно для групп).

## Чек-лист: что записать в отчёт
| № | Что записать | Где взять |
|---|---|---|
| 1 | Версии R, xcms, MSnbase | первая строка `module1_check/summary.csv`, `session_info_module1.txt` |
| 2 | Итог по каждому формату (CDF, mzML, mzXML): пройден / не пройден | `summary.csv` |
| 3 | Число признаков в таблице пиков по каждому формату (должно совпадать) | `summary.csv`, строка «одинаковое число признаков» |
| 4 | Время модуля 1 и пиковая память R | `summary.csv`; подробно — `project_faahKO_mzML/logs/modules_*.csv` |
| 5 | Строки журнала модуля 1: «Детектирование пиков», «Выравнивание», «Таблица пиков: N признаков, пропусков X %» | `project_faahKO_mzML/logs/run_*.log` |
| 6 | Что отчёт сформирован (модули 0–9 на mzML) | `project_faahKO_mzML/09_report/report.html` |
| 7 | Чтение mzML/mzXML из msdata | строки «msdata: …» в `summary.csv` |
| 8 | Ошибки и предупреждения, если были | консоль, `logs/run_*.log` (строки «ОШИБКА») |
| 9 | ОС, процессор, ОЗУ, время установки пакетов | вручную |

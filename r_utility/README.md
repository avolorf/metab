# UroMetabR 1.0.0 — утилита для полуавтоматической обработки данных ВЭЖХ-МС мочи

Контрольный экспериментальный образец программного инструмента (2 этап НИР). Язык R (≥ 4.0.0), работа локально,
без серверной инфраструктуры. Алгоритм совпадает с использованным в статье для журнала *Analytica*
(Python-реализация в `../pipeline/`); совпадение проверено на контрольных наборах (`tests/crosscheck.R`).

## Установка
1. Установить R ≥ 4.0.0 (и при желании RStudio).
2. Распаковать архив `UroMetabR_1.0.0.zip`, открыть каталог `r_utility`.
3. Один раз выполнить `Rscript install_packages.R` (или открыть файл в RStudio и нажать *Source*).

## Запуск
```bash
Rscript run_analysis.R config/config_example.yaml           # контрольный набор из example/
Rscript run_analysis.R config/config_example.yaml --quick   # быстрая проверка (3 повтора CV, 20 перестановок)
Rscript run_analysis.R config/config_s266.yaml              # батч s266 (нужен ../data/s266_table_for_calc.csv)
```
В RStudio: `setwd("…/r_utility"); source("run_analysis.R"); main("config/config_example.yaml")`.

## Входные данные
| файл | содержание |
|---|---|
| таблица признаков | `format: wide` — CSV/TSV, строки — инъекции (`Sample`), столбцы — признаки `RT_mz_аннотация` (или `mz / RT` из xcms); `format: msdial` — экспорт выравнивания MS-DIAL (Area/Height, txt); `format: xcms` — каталог mzML/mzXML/CDF (нужен Bioconductor-пакет xcms; модуль `R/01b_xcms.R`) |
| `metadata.csv` | `sample, group, sex (M/F или М/Ж), age, include` [, `exclude_reason`, `injection_order`] |
| ряд экзогенных соединений (необяз.) | `label, formula, ion, rt_min, in_index` → индекс загрязнения (ковариата) |
| ключевые слова экспозиций (необяз.) | `pattern, category` (drug / food / tobacco) |

QC и бланки распознаются по имени (`qc_pattern`, `blank_pattern`) или по строке *File type* экспорта MS-DIAL.

## Этапы обработки
1. Фильтры: детекция ≥ 80 % в QC и ≥ 50 % хотя бы в одной группе; QC/бланк ≥ 5 (если есть бланки); RSD в QC ≤ 30 %.
2. Дрейф: проверка leave-one-QC-out (медиана |Δlog2| без коррекции и с LOWESS); коррекция `none | lowess | auto`.
3. Нормализация PQN (референс — медианный спектр QC) и MSTUS для сравнения; импутация ½ минимума; log2.
4. Группировка ионов в соединения (ΔRT ≤ 0.02 мин, ρ Спирмена ≥ 0.8), подписи аддуктов/изотопов/потерь.
5. Индекс загрязнения (средний z-score членов ряда), категории экспозиций — исключаются из «находок» и моделей.
6. Для каждого сравнения: модель `log2 ~ группа + пол + возраст + индекс` (BH), подбор пар по возрасту внутри пола
   (оптимальное назначение, ±5 лет), сравнение с контролями с низким индексом; правило «находки».
7. Классификация: 20 × 5-блочная стратифицированная CV, отбор top-5 признаков **внутри** обучающего блока,
   логистическая регрессия (L2) и PLS-DA, базовые модели «только пол» / «только индекс», перестановочный тест.
8. Рисунки, таблицы (CSV + `results.xlsx`), HTML-отчёт, журнал, `session_info.txt`, контрольные суммы `manifest_md5.csv`.

## Результаты (`output`)
`report.html` — отчёт с рисунками; `results.xlsx` и `tables/*.csv` — все таблицы (журнал, демография, сводка,
находки, классификация, статистика по каждому сравнению, признаки с m/z и RT — база признаков, матрица log2);
`figures/*.png` — QC, распределения по пробам, PCA, индекс загрязнения, вулканы, тепловые карты, диаграммы
размаха, ROC, PQN vs MSTUS.

## Проверка корректности
```bash
Rscript tests/run_tests.R            # 19 тестов / 78 проверок: эталоны R/CRAN + контрольные наборы с известной истиной
Rscript tests/validation_metrics.R 5 # метрики на 5 парах контрольных наборов → tests/validation_metrics.csv
Rscript tests/crosscheck.R           # сверка с Python-реализацией (нужен python3 + numpy, pandas, scipy, statsmodels)
Rscript tests/make_control_data.R example/control_set_2 2   # новый контрольный набор (seed 2)
```

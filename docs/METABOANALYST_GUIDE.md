# Работа с MetaboAnalyst и базами данных: пошаговая инструкция

Все файлы лежат в `results/external/` (их создаёт `annotate_visualize.py`). Номер в начале имени файла совпадает с номером шага ниже.

Сайт: https://www.metaboanalyst.ca → **Click here to start**. Названия пунктов меню соответствуют версии 6.0; в других версиях формулировки могут немного отличаться.

## Важно перед началом
* Данные в файлах **уже обработаны**: коррекция дрейфа по QC, фильтр RSD_QC ≤ 30 %, нормализация PQN, импутация пропусков. Поэтому в MetaboAnalyst **не нужно** повторно нормализовать образцы (Sample normalization = *None*). Нужно только log-преобразование и шкалирование.
* Названия признаков имеют вид `m/z__RT(мин)`, например `583.2553__13.476`. Метки групп — `BC`, `PC`, `Control`.
* Для каждого анализа есть два варианта таблицы:
  * `_all7916` — все признаки, прошедшие фильтры качества;
  * `_IQRtop5000` — 5000 самых вариабельных признаков. Отбор по межквартильному размаху, без учёта групп, то есть так же, как фильтр IQR в самом MetaboAnalyst. **Используйте этот вариант:** MetaboAnalyst ограничивает число признаков и иначе потребует фильтрацию сам.
* Помните о двух систематических смещениях:
  1. **пол:** в группе рака в основном мужчины, в контроле в основном женщины;
  2. **источник проб:** во всех пробах «ДР» (в том числе почти во всех PC) есть ряд ПАВ-контаминантов.

  Шаг 2 позволяет учесть пол. Источник для PC учесть нельзя, потому что он почти полностью совпадает с диагнозом.

---

## Шаг 1. Statistical Analysis [one factor]: обзор данных, PCA, PLS-DA, тепловые карты
**Зачем:** увидеть, разделяются ли группы, найти отличающиеся признаки и построить графики для статьи.

1. На стартовой странице выберите модуль **Statistical Analysis [one factor]**.
2. Задайте параметры загрузки:
   * Data Type: **Peak intensities**;
   * Format: **Samples in rows (unpaired)**;
   * файл: `1_stats_BC_vs_Control_IQRtop5000.csv`.
   Нажмите *Submit*.
3. **Data Integrity Check.** Должно быть 74 образца, 2 группы (BC 29, Control 45), 5000 признаков, пропусков нет. Нажмите *Proceed*.
4. **Missing values.** Пропусков нет, оставьте по умолчанию.
5. **Data Filtering.** Выберите *None*: фильтрация уже сделана. Если без фильтра не пускает, выберите IQR.
6. **Normalization.** Параметры:
   * Sample normalization: **None**;
   * Data transformation: **Log transformation (base 10)**;
   * Data scaling: **Pareto scaling**.

   На графиках «до/после» распределения должны стать примерно симметричными.
7. Анализы (меню слева):
   * **Volcano plot**: FC threshold 2.0, t-test, p-value 0.05 с поправкой FDR. Это аналог нашего `figures/05_volcano_*`.
   * **PCA**: Scores (2D/3D), Loadings, Biplot. Сравните с `figures/04_pca_samples.png`.
   * **PLS-DA**: смотрите вкладки *Cross Validation* (Q2) и *Permutation* (поставьте 1000 перестановок). В статье указывайте Q2 и p-значение перестановочного теста, а не только график scores.
   * **OPLS-DA** (только для 2 групп): S-plot показывает признаки, которые сильнее всего отличают группы; дополнительно проверьте Permutation.
   * **Heatmap** (Cluster Analysis → Heatmaps):
     * Distance *Euclidean*, Clustering *Ward*, Standardization *Autoscale features*;
     * для читаемости включите *Use top 50 features by t-test/ANOVA*;
     * для карты средних по группам — *Display group averages*.
   * **Random Forest**: OOB error и важность признаков (*Variable importance*).
   * **Correlation / PatternHunter**: какие признаки коррелируют, например с биливердином.
8. Повторите шаги для `1_stats_PC_vs_Control_IQRtop5000.csv`.
9. Затем загрузите **трёхгрупповой** файл `1_stats_BC_PC_Control_IQRtop5000.csv`. В нём вместо t-теста **ANOVA** с post-hoc (Fisher LSD / Tukey), PCA и PLS-DA для трёх групп и тепловая карта средних (*Display group averages*).

> Ограничение: результаты шага 1 не учитывают пол.

## Шаг 2. Statistical Analysis [metadata table]: поправка на пол, влияние источника проб
**Зачем:** проверить, какие различия сохраняются после учёта пола, и увидеть эффект источника проб.

1. Выберите модуль **Statistical Analysis [metadata table]**.
2. Загрузите два файла:
   * данные: `2_metadata_module_data_IQRtop5000.csv` (образцы в строках; два контроля без указанного пола исключены, осталось 91 образец);
   * метаданные: `2_metadata_module_metadata.csv` — столбцы `Class` (BC/PC/Control), `Sex` (M/F) и `Source`. Значения `Source` — префиксы шифров в латинице: DR, OD, MNN, RR, RMP, KG, num.
3. На странице проверки метаданных отметьте все три столбца как **categorical**.
4. Normalization: None, Log10, Pareto (как в шаге 1).
5. Анализы:
   * **Linear Models with Covariate Adjustment** (limma):
     * Primary metadata: **Class**, reference group: **Control**;
     * Covariates: **Sex**;
     * сравнения *BC vs Control* и *PC vs Control*.

     Признаки, которые остаются значимыми (FDR < 0.05), не объясняются различием по полу. Сравните со столбцом `q_sex_adj` в нашей таблице.
   * **Interactive Heatmap / Heatmap with metadata**: включите полосы Class + Sex + Source. Будет видно, что ряд ПАВ (m/z 274.274, 302.306, 318.300…) окрашивает пробы по источнику DR, а не по диагнозу.
   * **iPCA / PCA с раскраской по Source**: проверка влияния источника проб.
   * **Не добавляйте `Source` как ковариату в сравнении PC vs Control.** Почти все PC — это DR, и модель станет неидентифицируемой. Именно поэтому PC в этом батче нельзя корректно отделить от эффекта источника.

## Шаг 3. Biomarker Analysis: корректная ROC
**Зачем:** честно оценить диагностическую способность, без завышения AUC.

1. Выберите модуль **Biomarker Analysis**. Загрузите `1_stats_BC_vs_Control_IQRtop5000.csv` с теми же параметрами, что в шаге 1. Нормализация: None / Log10 / Pareto (или Auto).
2. **Classical univariate ROC curve analyses** — AUC и 95 % ДИ каждого признака. Это описательная оценка. Если выбрать «лучшие» признаки по этой таблице и на тех же данных показывать их AUC, оценка будет завышена.
3. **Multivariate ROC curve exploration** — основной корректный анализ:
   * Classification method: **Linear SVM**, **PLS-DA** или **Random Forest** (сравните);
   * Feature ranking method: тот же метод или *Univariate t-statistic*;
   * MetaboAnalyst строит модели с 5, 10, 15, 25, 50, 100 признаками и **отбирает признаки внутри каждого цикла Monte-Carlo CV** (обучение на 2/3, тест на 1/3), поэтому AUC не завышен;
   * смотрите AUC с 95 % ДИ, *Predicted class probabilities*, *Feature importance* (как часто признак отбирается);
   * дополнительно запустите *Permutation test*.
4. **ROC curve based model creation & evaluation (Tester)** — ручная панель, например биливердин + LPC 16:0 + стероид C19H28O3. **Корректно только с hold-out:** отметьте часть образцов (~1/3) как *hold-out* **до** выбора признаков. Если панель уже выбрана по всем данным, AUC будет оптимистичным.
5. Повторите для `1_stats_PC_vs_Control_IQRtop5000.csv`.
6. Сверьте результаты с `results/roc_validation.xlsx`: там такая же схема (отбор признаков внутри CV), сделанная локально.

## Шаг 4. Functional Analysis — MS Peaks to Pathways (mummichog / GSEA)
**Зачем:** найти метаболические пути по m/z без точной идентификации соединений. Этот анализ мы уже сделали локально (`results/mummichog/`). MetaboAnalyst дополнительно даёт GSEA, библиотеку KEGG и интерактивную сеть.

1. Выберите модуль **Functional Analysis** (*MS Peaks to Pathways*).
2. Input: **Peak list profile** (m/z features with statistics). Файл: `3_peaks_BC_vs_Control.txt`, столбцы `m.z`, `p.value`, `t.score`, `r.t`. Пояснения:
   * p и знак эффекта взяты из модели с поправкой на пол;
   * ряд ПАВ исключён;
   * `t.score` = знак эффекта × −log10 p и служит только для ранжирования;
   * `r.t` указано в **секундах**.
3. Параметры:
   * Ion mode **Positive**, Mass tolerance **5 ppm**;
   * Retention time: **есть**, единицы **seconds**;
   * Rank features by: **P values**;
   * **Enforce primary ions**: да.
4. На следующей странице:
   * Algorithm: **mummichog** и **GSEA** (будет интегрированный результат);
   * p-value cutoff: **top 10 %** (по умолчанию);
   * Pathway library: **Homo sapiens (human) [MFN]** — та же модель, что у нас. Затем повторите с **Homo sapiens [KEGG]**.
5. Значимые пути: combined p < 0.05 и не меньше 3 найденных соединений. Если путь значим и в mummichog, и в GSEA, это надёжнее. Учтите:
   * p-значения mummichog не скорректированы на множественные сравнения;
   * большие пути (Tyrosine metabolism, ~170 соединений) часто дают ложноположительный результат из-за изомеров.
6. Повторите для `3_peaks_PC_vs_Control.txt` и `3_peaks_BC_vs_PC.txt`.

Вариант: в этот же модуль можно загрузить таблицу интенсивностей (`1_stats_*_IQRtop5000.csv`, Input = *Peak intensity table*). Тогда MetaboAnalyst сам посчитает t-тест, но без поправки на пол.

## Шаг 5. Enrichment Analysis (MSEA) по списку идентифицированных соединений
**Зачем:** проверить, обогащены ли курированные соединения известными путями или болезненными сигнатурами мочи.

1. Выберите модуль **Enrichment Analysis**, Input: **Compound list** (Over-representation analysis).
2. Вставьте содержимое `4_compounds_BC_HMDB.txt`, ID type: **HMDB ID**. Другой вариант — `4_compounds_BC_names.txt` с ID type **Compound name**.
3. Проверьте страницу **Name mapping**: все соединения должны сопоставиться.
4. Библиотеки:
   * **SMPDB** (метаболические пути);
   * **Disease signatures — Urine** (сравнение с известными сигнатурами заболеваний в моче);
   * **Drug pathways** — для лекарственных метаболитов.
5. Ограничения:
   * в списке BC около 15 соединений, это мало, поэтому результат поисковый;
   * без собственного референсного метаболома (*Reference metabolome*) фоном служат все соединения библиотеки, и обогащение получается немного завышенным.
6. Повторите для PC (`4_compounds_PC_*`). Там всего 3–4 соединения, поэтому для PC этот шаг практически неинформативен.

## Шаг 6. Pathway Analysis (KEGG, топология)
1. Выберите модуль **Pathway Analysis**, Input: **Compound list**, файл `4_compounds_BC_HMDB.txt`, ID type **HMDB**.
2. Параметры:
   * Pathway library: **KEGG — Homo sapiens**;
   * Enrichment method: **Hypergeometric test**;
   * Topology: **Relative-betweenness centrality**.
3. Результат — график Pathway Impact: по оси X влияние (impact), по оси Y −log10 p. Ожидаемо отклик дадут пути стероидных гормонов, глицерофосфолипидов, порфиринов (биливердин) и первичных желчных кислот.

## Шаг 7. KEGG Mapper — раскраска карт путей
1. Откройте https://www.kegg.jp/kegg/mapper/color.html.
2. Вставьте содержимое `6_kegg_mapper_BC.txt`: красный — повышено при BC, синий — понижено.
3. Search against: **hsa** (Homo sapiens) или **map**, затем *Exec*. Откройте карты, например *Porphyrin metabolism* (биливердин) и *Steroid hormone biosynthesis* (тетрагидрокортизон).

## Шаг 8. HMDB и METLIN — поиск кандидатов по массе
1. **HMDB LC-MS search** (https://hmdb.ca/spectra/ms/search):
   * вставьте `7_metlin_hmdb_mz.txt`;
   * Ion mode **Positive**;
   * Adducts **M+H, M+Na, M+NH4, M+H-H2O**;
   * Tolerance **5 ppm**.

   Если сервис ограничивает число масс, разбейте список на части.
2. **METLIN Batch Metabolite Search** (https://metlin.scripps.edu, нужна бесплатная регистрация): те же m/z, positive mode, аддукты [M+H]+, [M+Na]+, [M+NH4]+, 5 ppm.
3. Совпадение по массе — это только **кандидат** (MSI уровень 3). Для идентификации нужно:
   * извлечь из исходных данных MS/MS-спектры этих признаков и сравнить их со спектрами METLIN MS/MS, mzCloud или MoNA/GNPS (уровень 2);
   * подтвердить ключевые соединения стандартами по RT и MS/MS (уровень 1). В первую очередь это биливердин, стероид C19H28O3, тетрагидрокортизон и LPC.

## Рекомендуемый порядок
1. Шаг 1 (обзор) и шаг 2 (поправка на пол): какие различия реальны.
2. Шаг 3: насколько хорошо они разделяют группы (корректная ROC).
3. Шаги 4–7: биологическая интерпретация.
4. Шаг 8 и MS/MS: подтверждение идентичности.

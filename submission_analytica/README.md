# Пакет для подачи в *Analytica* (MDPI)
«Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study»

Собран из тех же расчётов (`results/`), что полная версия (`results/publication/`) и v4 (`publication_v4/`).
Пересобрать рисунки и таблицы: `python make_figures.py --variant v5` (папка `0_manuscript/` не трогается).

| папка / файл | содержимое |
|---|---|
| `0_manuscript/Methods_Results_v5.3_submission.md` | Methods + Results and Discussion + Conclusions (v5.3 с исправленными номерами таблиц и числами после итогового пересчёта) |
| `1_main_figures/` | Figures 1–4: TIFF (LZW) и PNG 600 dpi, PDF (вектор) |
| `2_main_tables/` | `Tables_main.xlsx` (Table 1 + сноска, Table 1 tests, Table 2 + сноска об изомерах), CSV-копии |
| `3_supplementary/` | Figures S1–S6; `Supplementary_Tables.xlsx` (S1–S7, см. ниже); CSV-копии |
| `figure_captions.md`, `table_notes.md` | подписи и сноски (англ.) |

Supplementary: Figure S1 тепловая карта; S2 возраст и пары; S3 PQN vs MSTUS; S4 пути (mummichog); S5 бланки
растворителя; S6 гомологический ряд и индекс. Table S1 все находки РМП; S2 классификация (с примечанием о РПЖ vs
контроль); S3 курированные аннотации (S3b — несходящаяся масса); S4 «есть/нет»; S5 статистика по соединениям
(4 листа); S6 анализы чувствительности; S7 проверка бланков.

## Что осталось сделать авторам
- [ ] Заполнить все места в [квадратных скобках] в тексте (сбор мочи, этика, LC–HRMS параметры). Параметры MS-DIAL вписаны в 2.3 по Таблице 3 отчёта НКР — подтвердить, что они относятся к батчу s266, и сверить английские названия метрик сравнения (simple/weighted dot product) с окном MS-DIAL.
- [ ] **Прибор:** в тексте указан Orbitrap Exploris 120, а в статье ABC — Orbitrap Fusion Lumos Tribrid; привести к фактическому прибору этого батча.
- [ ] В тексте сослаться на Table S3 и Table S5 (MDPI требует ссылки на все материалы приложения).
- [ ] Introduction, Abstract, ключевые слова, Author Contributions, Funding, IRB/Informed Consent Statements, Data Availability, Conflicts of Interest — по шаблону MDPI.
- [ ] Собрать рукопись в шаблоне MDPI (Word/LaTeX); таблицы 1–2 перенести в текст; Supplementary — один PDF с рисунками S1–S6 + `Supplementary_Tables.xlsx`.
- [ ] Сверить требования к рисункам в Instructions for Authors журнала (формат, разрешение, ширина).

## Обновление 15.10.2026
Внесены уточнённые пол и возраст (№20 — М, 68 лет; №62 — М, 59; №94 — Ж, 18) и выполнен полный пересчёт.
Актуальны только числа в `0_manuscript/` и в этой папке; тексты в `publication_v4/text/` и `docs/` относятся
к расчёту до уточнения.

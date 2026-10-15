# Вариант v4 — «Urinary Metabolomic Profiles of Bladder and Prostate Cancer by LC–HRMS: A Comparative Study»

Отдельная папка, чтобы можно было сравнить с полной версией (`results/publication/`, `docs/Results_and_Discussion.md`).
Расчёты те же (`results/`), меняются только оформление, терминология и состав рисунков.

* `text/Results_and_Discussion_v4.md` — раздел Results and Discussion (v4).
* `figures_tables/` — рисунки (PDF, TIFF, PNG 600 dpi), `Tables_main.xlsx`, `Supplementary_Tables.xlsx`, CSV-копии,
  `figure_captions.md`, `table_notes.md`. Пересобрать: `python make_figures.py --variant v4` (папка `text/` не трогается).

| v4 | что это | в полной версии |
|---|---|---|
| Figure 1 (a–c) | качество данных, без панели d | Figure 1 |
| Figure 2 | вулканы | Figure 3 |
| Figure 3 | 7 аннотированных находок | Figure 4 |
| Figure 4 | ROC: все пробы, с поправкой на индекс, пары по возрасту, базовая линия «пол» | Figure 5 |
| Figure 5 (необязательно) | маркеры предыдущего исследования | Figure 6 |
| Figure S1–S4 | тепловая карта, возраст/пары, PQN vs MSTUS, пути | те же |
| Figure S5 | бланки растворителя | Figure S5 |
| Figure S6 | гомологический ряд и индекс по группам | Figure 2 |
| Tables 1–2, S1–S8 | те же; «contamination index» → «exogenous-signal index», «stage 1» → «previous study» | те же |

Если сравнение с предыдущим исследованием не включается: убрать раздел 3.6, Figure 5 и Table S3 (лист S3_previous_markers).

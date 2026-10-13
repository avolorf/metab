"""Все параметры обработки в одном месте. Изменение этого файла или metadata/metadata.csv
и повторный запуск `python run_all.py` полностью воспроизводят результаты."""

# ---------------------------------------------------------------- пути
DATA = "data/s266_table_for_calc.csv"          # таблица признаков (Sample, Group, Diag, RT_mz_annotation...)
METADATA = "metadata/metadata.csv"             # sample, group, sex, age, include, exclude_reason, comment
OUT = "results"
CURATION = "curation/candidates_curated.csv"   # ручная курация кандидатов (формула, ион, ID, категория, MS2)
CONTAMINANT_SERIES = "curation/contaminant_series.csv"
EXPOSURE_KEYWORDS = "curation/exposure_keywords.csv"

SEED = 20261013
N_JOBS = -1                                    # параллельные перестановки (joblib)

# ---------------------------------------------------------------- дизайн
STUDY_GROUPS = ["Control", "BC", "PC"]
GROUP_LABELS = {"Control": "Контроль", "BC": "РМП", "PC": "РПЖ"}
QC_PREFIX = "QC"
BLANK_REGEX = r"^(blank|бланк|blk)"           # пробы, чьё имя или группа совпадает, считаются бланками
QC_EVERY = 20                                  # QC1 перед образцом 1, далее QC после каждых 20 инъекций

# сравнения: case / reference, ограничение по полу (None = оба пола)
COMPARISONS = {
    "BC_vs_Control": dict(case="BC", ref="Control", sex=None, primary=True),
    "BC_vs_PC": dict(case="BC", ref="PC", sex=None, primary=True),
    "PC_vs_Control": dict(case="PC", ref="Control", sex=None, primary=False),
    "PC_vs_maleControl": dict(case="PC", ref="Control", sex="M", primary=False),
}
COVARIATES = ["sex", "age", "contam_index"]    # основная модель: log2 ~ группа + ковариаты (complete case)

# подбор по возрасту: 1:1, ближайший возраст, допуск ±CALIPER лет, внутри пола (оптимальное назначение)
MATCH_CALIPER_YEARS = 5
MATCH_WITHIN_SEX = True

# ---------------------------------------------------------------- предобработка
MIN_QC_DETECT = 0.8
MIN_GROUP_DETECT = 0.5
MAX_QC_RSD = 0.30
BLANK_RATIO = 5.0                              # медиана QC / медиана бланков >= 5
DRIFT = "none"                                 # "none" | "lowess" (6 QC: коррекция не даёт выигрыша)
NORMALIZATION = "PQN"                          # основной метод: "PQN" | "MSTUS"
COMPARE_NORMALIZATIONS = ["PQN", "MSTUS"]
LOW_DILUTION_FACTOR = 0.2                      # PQN-фактор ниже -> проба помечается для проверки
PRESENCE_ABSENCE = (0.5, 0.8)                  # «есть/нет»: в одной группе детекция < 0.5, в другой >= 0.8

# ---------------------------------------------------------------- группировка признаков в соединения
GROUP_RT_TOL = 0.02                            # мин
GROUP_MIN_CORR = 0.8                           # Спирмен по образцам
PPM = 5.0

# ---------------------------------------------------------------- индекс загрязнения
CONTAM_SPLIT = "median"                        # высокий/низкий индекс: по медиане включённых проб

# ---------------------------------------------------------------- статистика и правила «находки»
FDR = 0.05
CONFIRM_P = 0.05                               # подтверждающие проверки (подбор по возрасту, низкий индекс)
CTRL_INDEX_P = 0.05                            # контроли high vs low индекс: p >= -> нет связи с загрязнением
CV_REPEATS = 20
CV_FOLDS = 5
PANEL_K = 5
N_PERM = 1000
N_BOOT = 2000
EXCLUDE_DRUG_LIKE_UNKNOWN_FROM_PANEL = True    # низкоуверенные «лекарственные» аннотации -> не в панель
PANEL_MAX_ABS_RHO_CONTAM = 0.5                 # соединения с |ρ(индекс загрязнения)| ≥ порога -> не в панель
MIN_PAIRS_CV = 10                              # CV по подобранным парам только при ≥ 10 парах (иначе «не оценивается»)

# ---------------------------------------------------------------- сверка с маркерами первого этапа / статьи ABC
PRIOR_MARKERS = "curation/prior_markers.csv"   # если файла нет — раздел пропускается
PRIOR_RT_TOL = 0.2                             # мин; при другом градиенте оставить rt_min пустым (только масса)

# ---------------------------------------------------------------- пути метаболизма (приложение)
RUN_PATHWAYS = True
MUMMICHOG_TOP_FRAC = 0.10
MUMMICHOG_PERM = 500

FIG_DPI = 300

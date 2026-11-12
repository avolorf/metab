#!/usr/bin/env python3
"""Утилита запуска обработки ВЭЖХ-МСВР данных (фреймворк поиска биомаркеров).

Без аргументов открывает меню. Команды:

    python metab.py check                      проверить Python, пакеты и наличие исходных файлов
    python metab.py convert ЭКСПОРТ.txt --out data/новая/таблица.csv
                                               экспорт выравнивания MS-DIAL -> таблица для фреймворка
    python metab.py run s107                   полная обработка выборки s107 + таблицы, рисунки, проверка ABC
    python metab.py run s266                   полная обработка выборки s266 + анализы чувствительности + рисунки статьи
    python metab.py run s107 --quick           быстрый пробный прогон (в results_quick_s107/, основные результаты не трогаются)
    python metab.py chapter                    пересобрать главу диссертации (docx) по готовым результатам
    python metab.py all                        s266 -> s107 -> глава (полный цикл)
    python metab.py new s300 --export ЭКСПОРТ.txt
                                               новая выборка: таблица, заготовка метаданных, настройки datasets/s300.json
    python metab.py run s300                   обработка новой выборки (после заполнения metadata/metadata_s300.csv)

Журналы каждого шага пишутся в logs/ (имя файла = дата_время_шаг.log).
"""
import argparse
import importlib.metadata as md
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LOGS = ROOT / "logs"

DATASETS = {
    "s266": {
        "title": "s266: РМП, РПЖ и контроль (94 участника)",
        "required": ["data/s266_table_for_calc.csv", "metadata/metadata_stage2.csv", "curation/candidates_curated.csv"],
        "optional": {"data/raw_check": "mzML для проверки бланков (blank_check.py) — без них шаг пропускается"},
        "out": "results",
    },
    "s107": {
        "title": "s107: РМП и контроль (107 проб, в анализе 104)",
        "required": ["data/s107/s107_table.csv", "metadata/metadata_s107.csv", "data/s107/table_nastya.csv",
                     "data/s107/res_annotated_2.csv", "curation_s107/candidates_curated.csv",
                     "results/tables/stats_BC_vs_Control.csv"],
        "optional": {"data/abc/T1_ABC_markers_in_s107.xlsx": "таблица соединений статьи ABC — без неё шаг проверки ABC пропускается"},
        "out": "results_s107",
    },
}


def say(msg=""):
    print(msg, flush=True)


def step(name, cmd, env_extra=None):
    """Запустить шаг, показать вывод на экране и сохранить его в logs/."""
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{datetime.now():%Y-%m-%d_%H-%M-%S}_{name}.log"
    env = dict(os.environ, PYTHONIOENCODING="utf-8", MPLBACKEND="Agg", **(env_extra or {}))
    say(f"\n=== {name}: {' '.join(str(c) for c in cmd)}")
    say(f"    журнал: {log.relative_to(ROOT)}")
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.Popen([sys.executable, *map(str, cmd)], cwd=ROOT, env=env, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        for line in p.stdout:
            fh.write(line)
            if not any(w in line for w in ("Warning", "warnings.warn", "findfont")):
                print("    " + line.rstrip(), flush=True)
        rc = p.wait()
    dt = time.time() - t0
    if rc != 0:
        say(f"!!! Шаг «{name}» завершился с ошибкой (код {rc}) через {dt / 60:.1f} мин. Подробности — в {log.name}")
        sys.exit(rc)
    say(f"=== {name}: готово за {dt / 60:.1f} мин")


# ------------------------------------------------------------------ check
def check(verbose=True):
    ok = True
    v = sys.version_info
    say(f"Python {v.major}.{v.minor}.{v.micro} ({sys.executable})")
    if v < (3, 11):
        say("  !! нужен Python 3.11 или новее (рекомендуется 3.13)")
        ok = False
    say("\nПакеты (requirements.txt):")
    ready = []
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, _, want = line.partition("==")
        try:
            have = md.version(name)
            mark = "ok" if have == want else f"установлена {have} (в requirements {want}) — обычно работает"
        except md.PackageNotFoundError:
            mark, ok = "НЕ УСТАНОВЛЕН -> pip install -r requirements.txt", False
        say(f"  {name:14s} {mark}")
    say("\nИсходные файлы:")
    for ds, d in DATASETS.items():
        miss = [f for f in d["required"] if not (ROOT / f).exists()]
        say(f"  {d['title']}: {'все на месте' if not miss else 'НЕ ХВАТАЕТ ' + ', '.join(miss)}")
        if not miss:
            ready.append(ds)
        for f, why in d["optional"].items():
            if not (ROOT / f).exists():
                say(f"     (нет {f}: {why})")
    try:
        subprocess.run(["pandoc", "--version"], capture_output=True, check=True)
        say("\npandoc: найден (нужен только для сборки главы в docx)")
    except (OSError, subprocess.CalledProcessError):
        say("\npandoc: не найден — глава соберётся только в .md (для docx установите pandoc: https://pandoc.org)")
    if not ok:
        say("\nИтог: сначала установите недостающее (см. выше)")
    elif ready:
        say("\nИтог: установка в порядке; к обработке готовы: " + ", ".join(ready))
    else:
        say("\nИтог: установка в порядке, но исходных данных нет — скопируйте их в data/ (docs/Запуск_на_ПК.md, раздел 3)")
    return ok


def need(ds):
    miss = [f for f in DATASETS[ds]["required"] if not (ROOT / f).exists()]
    if miss:
        say("Не хватает файлов: " + ", ".join(miss) + "\n(исходные данные не хранятся в git — скопируйте их в папку data/, "
            "см. docs/Запуск_на_ПК.md, раздел 3)")
        sys.exit(1)


# ------------------------------------------------------------------ run
def run(ds, quick=False):
    if ds not in DATASETS:
        js = ROOT / "datasets" / f"{ds}.json"
        if not js.exists():
            say(f"Неизвестная выборка {ds}. Известные: {', '.join(list(DATASETS) + custom_datasets())}")
            sys.exit(1)
        import json
        st = json.loads(js.read_text(encoding="utf-8"))
        for f in (st["DATA"], st["METADATA"]):
            if not (ROOT / f).exists():
                say(f"Нет файла {f}")
                sys.exit(1)
        env = {"METAB_DATASET": ds, "METAB_SETTINGS": str(js)}
        out = f"results_quick_{ds}" if quick else st["OUT"]
        step(f"{ds}_run_all", ["run_all.py", "--out", out] + (["--quick"] if quick else []), env)
        say(f"\nВыборка {ds} обработана. Результаты: {out}/ (results.xlsx, tables/, figures/, log.txt)")
        return
    need(ds)
    env = {"METAB_DATASET": ds}
    if quick:
        out = f"results_quick_{ds}"
        say(f"Быстрый пробный прогон {ds}: 20 перестановок, 3 повтора CV, без путей; результаты в {out}/ "
            "(основные результаты не изменяются)")
        step(f"{ds}_quick", ["run_all.py", "--quick", "--out", out], env)
        return
    t0 = time.time()
    step(f"{ds}_1_run_all", ["run_all.py"], env)
    if ds == "s266":
        if (ROOT / "data/raw_check").exists():
            step("s266_2_blank_check", ["blank_check.py"], env)
        step("s266_3_sensitivity_imputation", ["sensitivity_imputation.py"], env)
        if (ROOT / "results/blank_check").exists():
            step("s266_4_sensitivity_index", ["sensitivity_index.py"], env)
        step("s266_5_figures_submission", ["make_figures.py", "--variant", "v5"], env)
    else:
        step("s107_2_report", ["report_s107.py"], env)
        if (ROOT / "data/abc/T1_ABC_markers_in_s107.xlsx").exists():
            step("s107_3_abc_table5", ["report_s107/abc_table5_check.py"], env)
    say(f"\nВыборка {ds} обработана за {(time.time() - t0) / 60:.0f} мин. Результаты: {DATASETS[ds]['out']}/"
        + (" и report_s107/" if ds == "s107" else " и submission_analytica/"))


def chapter():
    for f in ["results/tables/stats_BC_vs_Control.csv", "results_s107/tables/stats_BC_vs_Control.csv",
              "report_s107/Tables_s107.xlsx"]:
        if not (ROOT / f).exists():
            say(f"Нет {f}: сначала выполните обработку обеих выборок")
            sys.exit(1)
    step("chapter", ["report_chapter_ru/build_chapter.py"])


def convert(export, out, pattern=None):
    cmd = ["tools/msdial_to_table.py", export, "--out", out]
    if pattern:
        cmd += ["--sample-pattern", pattern]
    step("convert", cmd)


# ------------------------------------------------------------------ новая выборка
def custom_datasets():
    d = ROOT / "datasets"
    return sorted(p.stem for p in d.glob("*.json")) if d.exists() else []


CONTROL_WORDS = ("control", "контрол", "healthy", "здоров", "norm", "кг")


def new(name, export=None, table=None, control=None, sample_pattern=None, qc_every=None, batches=None):
    """Подготовить новую выборку: таблица, метаданные, курация, настройки datasets/ИМЯ.json."""
    import json
    import re
    import shutil
    import pandas as pd
    sys.path.insert(0, str(ROOT))
    if not re.fullmatch(r"[A-Za-z0-9_]+", name) or name in DATASETS:
        say("Имя выборки — латиница, цифры и «_», не s266/s107 (например s300)")
        sys.exit(1)
    ddir = ROOT / "data" / name
    ddir.mkdir(parents=True, exist_ok=True)
    tab = ddir / f"{name}_table.csv"
    if export:
        convert(export, str(tab.relative_to(ROOT)), sample_pattern)
    elif table:
        shutil.copy(table, tab)
    if not tab.exists():
        say(f"Нет таблицы {tab.relative_to(ROOT)}: укажите --export (экспорт MS-DIAL) или --table")
        sys.exit(1)
    stem = tab.with_suffix("")
    inj_p = Path(f"{stem}_injections.csv")
    inj = pd.read_csv(inj_p, dtype=str) if inj_p.exists() else None

    # метаданные: заготовка из конвертера, контрольная группа -> Control
    md_p = ROOT / "metadata" / f"metadata_{name}.csv"
    tpl_p = Path(f"{stem}_metadata_template.csv")
    if md_p.exists():
        md = pd.read_csv(md_p, dtype={"sample": str})
        say(f"metadata/{md_p.name} уже есть — не перезаписываю (группы беру из него)")
    else:
        if tpl_p.exists():
            md = pd.read_csv(tpl_p, dtype={"sample": str})
        else:
            t = pd.read_csv(tab, usecols=["Sample", "Group"], dtype=str)
            t = t[~t.Group.str.lower().isin(["qc", "blank"])]
            md = pd.DataFrame({"sample": t.Sample, "group": t.Group, "sex": "", "age": "", "include": 1,
                               "exclude_reason": "", "comment": ""})
        cls = md["group"].fillna("").astype(str).str.strip()
        ctrl = control or next((c for c in cls.unique() if c and c.lower().startswith(CONTROL_WORDS)), None)
        if not ctrl:
            say("Не удалось определить контрольную группу: укажите --control ИМЯ_КЛАССА (как в MS-DIAL). "
                f"Найдены классы: {sorted(set(cls))}")
            sys.exit(1)
        md["group"] = [("Control" if c == ctrl else re.sub(r"\W+", "_", c).strip("_")) for c in cls]
        md.to_csv(md_p, index=False)
        say(f"metadata/{md_p.name}: заготовка создана (контроль = класс «{ctrl}»). ЗАПОЛНИТЕ пол (M/F), возраст, "
            "include (1/0) и exclude_reason")
    groups = [g for g in md.loc[md.include.astype(str) == "1", "group"].dropna().astype(str).unique() if g]
    if "Control" not in groups:
        say("В метаданных нет группы Control — переименуйте контрольную группу в Control")
        sys.exit(1)
    cases = sorted(g for g in groups if g != "Control")
    comps = {f"{g}_vs_Control": dict(case=g, ref="Control", sex=None, primary=True) for g in cases}
    comps.update({f"{a}_vs_{b}": dict(case=a, ref=b, sex=None, primary=True)
                  for i, a in enumerate(cases) for b in cases[i + 1:]})

    # курация: пустой список кандидатов, MS2-совпадения из MS-DIAL, экзогенный ряд с выравниванием RT
    cdir = ROOT / f"curation_{name}"
    cdir.mkdir(exist_ok=True)
    cur = cdir / "candidates_curated.csv"
    if not cur.exists():
        cols = pd.read_csv(ROOT / "curation_s107" / "candidates_curated.csv", nrows=0).columns
        pd.DataFrame(columns=cols).to_csv(cur, index=False)
    ms2_p = Path(f"{stem}_ms2.csv")
    if ms2_p.exists():
        m = pd.read_csv(ms2_p, dtype=str)
        if "MS/MS matched" in m:
            m[m["MS/MS matched"].str.lower() == "true"][["feature"]].to_csv(cdir / "ms2_identified.csv", index=False)
    from pipeline.features import ADDUCTS, formula_mass
    from pipeline.io import parse_features
    info = parse_features(list(pd.read_csv(tab, nrows=0).columns.drop(["Sample", "Group"], errors="ignore")))
    ser = pd.read_csv(ROOT / "curation" / "contaminant_series.csv")
    shifts = []
    for _, r in ser.iterrows():
        th = formula_mass(r.formula) + ADDUCTS[r.ion]
        d = info[((info.mz - th).abs() / th * 1e6 <= 5) & ((info.rt_min - r.rt_min).abs() <= 2.0)]
        if len(d):
            shifts.append(float((d.rt_min - r.rt_min).iloc[(d.rt_min - r.rt_min).abs().argmin()]))
    shift = float(pd.Series(shifts).median()) if shifts else 0.0
    n_ok = sum(abs(x - shift) <= 0.1 for x in shifts)
    ser2 = ser.copy()
    ser2["rt_min"] = (ser2.rt_min + shift).round(3)
    ser2["comment"] = ser2.comment.fillna("").astype(str) + f" RT сдвинуто на {shift:+.2f} мин (metab.py new)"
    ser2.to_csv(cdir / "contaminant_series.csv", index=False)
    use_idx = n_ok >= 2
    say(f"Экзогенный ряд: найдено {n_ok} членов (сдвиг RT {shift:+.2f} мин относительно s266)"
        + ("" if use_idx else " — меньше 2, индекс экзогенного сигнала использоваться не будет"))

    # серии и порядок инъекций
    settings = {"_comment": f"Настройки выборки {name}; создано metab.py new {datetime.now():%Y-%m-%d}. "
                            "Можно править вручную (ключи — параметры pipeline/config.py).",
                "DATA": str(tab.relative_to(ROOT)), "METADATA": str(md_p.relative_to(ROOT)), "OUT": f"results_{name}",
                "CURATION": str(cur.relative_to(ROOT)), "CONTAMINANT_SERIES": f"curation_{name}/contaminant_series.csv",
                "MS2_IDENTIFIED": f"curation_{name}/ms2_identified.csv", "PRIOR_MARKERS": f"curation_{name}/prior_markers.csv",
                "RT_ANCHORS": f"curation_{name}/rt_anchors.csv", "STUDY_GROUPS": ["Control"] + cases,
                "GROUP_LABELS": {"Control": "Контроль", **{g: g for g in cases}}, "COMPARISONS": comps,
                "CONTAM_SPLIT": "control_median", "USE_CONTAM_INDEX": use_idx}
    if inj is not None and inj["injection_order"].fillna("").str.strip().ne("").all():
        settings["INJECTIONS"] = str(inj_p.relative_to(ROOT))
        b = sorted(inj["batch"].dropna().astype(str).str.strip().replace("", pd.NA).dropna().unique())
        settings["BATCHES"] = b if len(b) > 1 else None
        settings["BATCH_CORRECTION"] = len(b) > 1
        say(f"Порядок инъекций и серии — из проекта MS-DIAL ({len(b) or 1} сер.)")
    else:
        if batches:
            settings["BATCHES"] = [[int(x) for x in part.split("-")] for part in batches.split(",")]
            settings["BATCH_CORRECTION"] = len(settings["BATCHES"]) > 1
        settings["QC_EVERY"] = qc_every or 6
        say("Порядок инъекций восстанавливается по номерам проб: QC в начале серии и после каждых "
            f"{settings['QC_EVERY']} проб" + (f"; серии {batches}" if batches else "; одна серия"))
    (ROOT / "datasets").mkdir(exist_ok=True)
    (ROOT / "datasets" / f"{name}.json").write_text(json.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8")
    say(f"\nГотово: datasets/{name}.json. Группы: Control + {', '.join(cases)}; сравнения: {', '.join(comps)}.")
    say(f"Дальше: 1) заполните metadata/{md_p.name}; 2) python metab.py run {name}  (или пункт меню)")


# ------------------------------------------------------------------ меню
def menu():
    items = [("Проверить установку и исходные файлы", lambda: check()),
             ("Обработать выборку s107 (около 40–60 мин)", lambda: confirm("results_s107/") and run("s107")),
             ("Обработать выборку s266 (около 1–1.5 ч)", lambda: confirm("results/ и submission_analytica/") and run("s266")),
             ("Быстрый пробный прогон s107 (несколько минут)", lambda: run("s107", quick=True)),
             ("Пересобрать главу диссертации", chapter),
             ("Полный цикл: s266 -> s107 -> глава", lambda: confirm("results/, results_s107/, report_s107/") and full()),
             ("Конвертировать экспорт MS-DIAL в таблицу", ask_convert),
             ("Новая выборка: подготовить из экспорта MS-DIAL", ask_new),
             ("Обработать новую выборку", ask_run_custom)]
    while True:
        say("\n" + "=" * 60 + "\n  Обработка ВЭЖХ-МСВР данных — меню\n" + "=" * 60)
        for i, (t, _) in enumerate(items, 1):
            say(f"  {i}. {t}")
        say("  0. Выход")
        c = input("Выберите пункт: ").strip()
        if c in ("0", "", "q"):
            return
        if c.isdigit() and 1 <= int(c) <= len(items):
            try:
                items[int(c) - 1][1]()
            except SystemExit:
                pass
        input("\nНажмите Enter, чтобы вернуться в меню…")


def confirm(what):
    a = input(f"Результаты в {what} будут пересозданы. Продолжить? [д/н]: ").strip().lower()
    return a in ("д", "да", "y", "yes")


def ask_convert():
    exp = input("Путь к файлу экспорта MS-DIAL (.txt): ").strip().strip('"')
    out = input("Куда сохранить таблицу (например data/s300/s300_table.csv): ").strip().strip('"')
    pat = input("Шаблон номера образца в имени файла (Enter — не использовать; пример (\\d+)$): ").strip() or None
    convert(exp, out, pat)


def ask_new():
    name = input("Имя новой выборки (латиницей, например s300): ").strip()
    exp = input("Путь к файлу экспорта MS-DIAL (.txt): ").strip().strip('"')
    ctrl = input("Имя класса контрольной группы в MS-DIAL (Enter — определить автоматически): ").strip() or None
    pat = input("Шаблон номера пробы в имени файла (Enter — не использовать; пример (\\d+)$): ").strip() or None
    q = input("QC вводился после каждых N проб — N (Enter — 6; нужно, только если в MS-DIAL не заполнен "
              "Analytical order): ").strip()
    bt = input("Серии по номерам проб, например 1-40,41-80 (Enter — взять из MS-DIAL или одна серия): ").strip() or None
    new(name, export=exp, control=ctrl, sample_pattern=pat, qc_every=int(q) if q.isdigit() else None, batches=bt)


def ask_run_custom():
    ds = custom_datasets()
    if not ds:
        say("Новых выборок нет — сначала пункт «Новая выборка»")
        return
    say("Выборки: " + ", ".join(ds))
    name = input("Какую обработать: ").strip()
    if confirm(f"results_{name}/"):
        run(name)


def full():
    run("s266")
    run("s107")
    chapter()


def main():
    os.chdir(ROOT)
    if len(sys.argv) == 1:
        return menu()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="проверить установку и исходные файлы")
    r = sub.add_parser("run", help="обработать выборку")
    r.add_argument("dataset", help="s266, s107 или имя новой выборки (datasets/ИМЯ.json)")
    r.add_argument("--quick", action="store_true", help="быстрый пробный прогон в отдельную папку")
    sub.add_parser("chapter", help="пересобрать главу")
    n = sub.add_parser("new", help="подготовить новую выборку из экспорта MS-DIAL")
    n.add_argument("name", help="имя выборки латиницей, например s300")
    n.add_argument("--export", help="файл экспорта выравнивания MS-DIAL (.txt)")
    n.add_argument("--table", help="или уже готовая таблица (.csv: Sample, Group, признаки)")
    n.add_argument("--control", help="имя класса контрольной группы в MS-DIAL (если не распознаётся само)")
    n.add_argument("--sample-pattern", default=None, help="номер пробы из имени файла, например (\\d+)$")
    n.add_argument("--qc-every", type=int, help="QC после каждых N проб (если в проекте нет порядка инъекций)")
    n.add_argument("--batches", help="серии по номерам проб, например 1-36,37-72 (если в проекте нет Batch ID)")
    sub.add_parser("all", help="s266 -> s107 -> глава")
    c = sub.add_parser("convert", help="экспорт MS-DIAL -> таблица для фреймворка")
    c.add_argument("export")
    c.add_argument("--out", required=True)
    c.add_argument("--sample-pattern", default=None)
    a = ap.parse_args()
    if a.cmd == "check":
        sys.exit(0 if check() else 1)
    if a.cmd == "run":
        run(a.dataset, a.quick)
    elif a.cmd == "chapter":
        chapter()
    elif a.cmd == "all":
        full()
    elif a.cmd == "new":
        new(a.name, a.export, a.table, a.control, a.sample_pattern, a.qc_every, a.batches)
    elif a.cmd == "convert":
        convert(a.export, a.out, a.sample_pattern)


if __name__ == "__main__":
    main()

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


# ------------------------------------------------------------------ меню
def menu():
    items = [("Проверить установку и исходные файлы", lambda: check()),
             ("Обработать выборку s107 (около 40–60 мин)", lambda: confirm("results_s107/") and run("s107")),
             ("Обработать выборку s266 (около 1–1.5 ч)", lambda: confirm("results/ и submission_analytica/") and run("s266")),
             ("Быстрый пробный прогон s107 (несколько минут)", lambda: run("s107", quick=True)),
             ("Пересобрать главу диссертации", chapter),
             ("Полный цикл: s266 -> s107 -> глава", lambda: confirm("results/, results_s107/, report_s107/") and full()),
             ("Конвертировать экспорт MS-DIAL в таблицу", ask_convert)]
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
    r.add_argument("dataset", choices=list(DATASETS))
    r.add_argument("--quick", action="store_true", help="быстрый пробный прогон в отдельную папку")
    sub.add_parser("chapter", help="пересобрать главу")
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
    elif a.cmd == "convert":
        convert(a.export, a.out, a.sample_pattern)


if __name__ == "__main__":
    main()

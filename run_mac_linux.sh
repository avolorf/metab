#!/usr/bin/env bash
# Запуск обработки ВЭЖХ-МСВР данных:  bash run_mac_linux.sh   (меню)  или  bash run_mac_linux.sh run s107
# При первом запуске создаётся окружение .venv и устанавливаются пакеты (нужен интернет, 3-10 мин).
set -e
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
if [ ! -x .venv/bin/python ]; then
    echo "Первый запуск: создаю окружение Python в папке .venv ..."
    PY=$(command -v python3.13 || command -v python3)
    "$PY" -m venv .venv
    .venv/bin/python -m pip install --upgrade pip
    .venv/bin/python -m pip install -r requirements.txt
fi
exec .venv/bin/python metab.py "$@"

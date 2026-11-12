@echo off
rem Запуск обработки ВЭЖХ-МСВР данных (двойной щелчок по файлу).
rem При первом запуске создаётся окружение .venv и устанавливаются пакеты (нужен интернет, 3-10 мин).
chcp 65001 >nul
cd /d "%~dp0"
set PYTHONIOENCODING=utf-8
if not exist ".venv\Scripts\python.exe" (
    echo Первый запуск: создаю окружение Python в папке .venv ...
    py -3 -m venv .venv 2>nul || python -m venv .venv
    if errorlevel 1 (
        echo Не найден Python. Установите Python 3.13 с python.org ^(отметьте "Add python.exe to PATH"^) и запустите снова.
        pause
        exit /b 1
    )
    ".venv\Scripts\python.exe" -m pip install --upgrade pip
    ".venv\Scripts\python.exe" -m pip install -r requirements.txt
    if errorlevel 1 (
        echo Не удалось установить пакеты. См. сообщения выше и раздел "Типичные проблемы" в docs\Запуск_на_ПК.md
        pause
        exit /b 1
    )
)
".venv\Scripts\python.exe" metab.py %*
pause

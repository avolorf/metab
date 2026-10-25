@echo off
chcp 65001 >nul
cd /d "%~dp0"
set "RS="
for /f "delims=" %%i in ('where Rscript 2^>nul') do if not defined RS set "RS=%%i"
if not defined RS for /d %%d in ("%ProgramFiles%\R\R-*") do if exist "%%d\bin\Rscript.exe" set "RS=%%d\bin\Rscript.exe"
if not defined RS (
  echo R не найден. Установите R версии 4.0.0 или новее: https://cran.r-project.org/bin/windows/base/
  pause
  exit /b 1
)
echo Запуск UroMetabR: "%RS%"
"%RS%" --vanilla app\launch.R
if errorlevel 1 pause

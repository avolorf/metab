#!/bin/sh
# Запуск графического интерфейса UroMetabR (Linux, macOS)
cd "$(dirname "$0")" && Rscript --vanilla app/launch.R

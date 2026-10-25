#!/bin/sh
# Выпуск версии: листинг для депонирования, архив исходного кода, документы. Запуск из папки утилиты:
#   sh make_release.sh            (результат — в ../dist)
set -e
V=$(sed -n 's/^UROMETABR_VERSION <- "\(.*\)"/\1/p' R/00_core.R)
D=$(sed -n 's/^UROMETABR_DATE <- "\(.*\)"/\1/p' R/00_core.R)
OUT=../dist; mkdir -p "$OUT"
FILES="run_analysis.R $(ls R/*.R) app/app.R app/launch.R templates/report.Rmd Start_UroMetabR.bat start_urometabr.sh Create_desktop_shortcut.vbs"
L="$OUT/UroMetabR_v${V}_listing.txt"
{ echo "Программный инструмент для неинвазивной вспомогательной диагностики рака мочеполовой системы на основе анализа мочи методом ВЭЖХ-МС"
  echo "UroMetabR, версия $V от $D. Листинг исходного кода."
  for f in $FILES; do echo; echo "==================== $f ===================="; tr -d '\r' < "$f" | awk '{printf "%5d  %s\n", NR, $0}'; done; } > "$L"
Z="$OUT/UroMetabR_v${V}_${D}.zip"; rm -f "$Z"
cd .. && zip -qr "dist/$(basename "$Z")" r_utility -x "r_utility/tests/ispytaniya_results/runs/*" "r_utility/docs/*.docx" "r_utility/make_release.sh" && cd r_utility
node docs/build_docx.js docs/content_komplekt.js "$OUT/UroMetabR_v${V}_komplekt_registracii.docx" "Комплект для фиксации программы"
node docs/build_docx.js docs/content_zapiska.js "$OUT/UroMetabR_v${V}_poyasnitelnaya_zapiska.docx" "UroMetabR: пояснительная записка и испытания"
ls -la "$OUT"

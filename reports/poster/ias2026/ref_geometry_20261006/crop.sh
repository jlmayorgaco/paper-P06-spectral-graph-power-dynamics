#!/bin/bash
# usage: crop.sh name x0 y0 x1 y1 (ref px) shiftY
n=$1; x0=$2; y0=$3; x1=$4; y1=$5; s=${6:-0}
X=$(python -c "print(int(4.651*$x0))"); Y=$(python -c "print(int(4.649*($y0+$s)))")
W=$(python -c "print(int(4.651*($x1-$x0)))"); H=$(python -c "print(int(4.649*($y1-$y0)))")
pdftoppm -r 150 -png -singlefile -x $X -y $Y -W $W -H $H ${PDF:-build_polish/main.pdf} ${OUT:-crops_polish}/$n 2>/dev/null

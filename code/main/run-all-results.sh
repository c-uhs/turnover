#!/bin/bash

# run this script from: turnover/code/

# ./main/sensitivity-run.sh
for todo in \
  phi-base         \
  fit              \
  compare-tpaf     \
  sensitivity-plot \
  flows            \
; do
  echo "##################################################"
  echo $todo
  python3 main/main.py $todo
done
echo "##################################################"

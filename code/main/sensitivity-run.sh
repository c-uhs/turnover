#!/bin/bash

# run this script from: turnover/code/

i1=${1:-0}
i2=${2:-31}
echo "=================================================="
for i in `seq $i1 $i2`; do
  echo "--------------------------------------------------"
  echo "submitting run #$i"
  echo """#!/bin/bash
# module load anaconda3
python3 main/main.py sensitivity-run $i
""" > main/jobs/job-$i.sh
  chmod +x main/jobs/job-$i.sh
  # parallel (scinet)
  # sbatch --nodes=1 --chdir=. --time=0:15:00 ./main/jobs/job-$i.sh
  # serial (local)
  ./main/jobs/job-$i.sh > ./main/jobs/job-$i.out &
  sleep 5
done
echo "=================================================="
echo "done"

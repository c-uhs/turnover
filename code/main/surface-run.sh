#!/bin/bash

# run this script from: swazi-fsw/code/models/main/fit/scinet/

# $1 (argument #1) = number of jobs
# seeds selected sequentially from 1
for ((i=0;i<15;i++)); do
  echo "=================================================="
  echo "submitting opt #$i"
  echo "=================================================="
  echo """#!/bin/bash
module load anaconda3
python main/main.py surface-run $i
""" > main/jobs/job-$i.sh
  chmod +x main/jobs/job-$i.sh
  # parallel (scinet)
  sbatch --nodes=1 --chdir=. --time=0:15:00 ./main/jobs/job-$i.sh
done
echo "=================================================="
echo "done"
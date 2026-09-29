#!/bin/bash
#SBATCH --job-name=symcheck
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=2G
#SBATCH --time=04:30:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/sym_check.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/sym_check.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Symmetric HJ(4,2) calibration at n=12,13 (known SAT): kissat vs incidence+avoid vs dynseek+avoid.
# 18 runs (n=12,13 x 3 arms x 3 seeds), 9 at a time: the 9 fast n=12 runs, then the 9 long n=13 runs.
# Measured ~30 MB per run, so 2G is ample.
# Set KISSAT to your kissat binary if it is not on PATH (the kissat arm is skipped if not found).
KISSAT=${KISSAT:-$(command -v kissat)}
cd /home/hrmiller/CDCL-Experimentation
python3 -u experiments/hj_sym_check.py --n 12,13 --seeds 3 --timeout 14400 --workers 9 \
  --kissat "${KISSAT:-kissat}" --out experiments/results/sym_check.csv --resume

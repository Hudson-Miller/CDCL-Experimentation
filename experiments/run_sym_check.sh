#!/bin/bash
#SBATCH --job-name=symcheck
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=04:30:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/sym_check.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/sym_check.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Symmetric HJ(4,2) calibration at n=12,13 (known SAT): kissat vs incidence+avoid vs dynseek+avoid.
# Set KISSAT to your kissat binary if it is not on PATH (the kissat arm is skipped if not found).
KISSAT=${KISSAT:-$(command -v kissat)}
cd /home/hrmiller/CDCL-Experimentation
python3 -u experiments/hj_sym_check.py --n 12,13 --seeds 3 --timeout 14400 --workers 24 \
  --kissat "${KISSAT:-kissat}" --out experiments/results/sym_check.csv --resume

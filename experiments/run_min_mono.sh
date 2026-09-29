#!/bin/bash
#SBATCH --job-name=minmono
#SBATCH --partition=cpu
#SBATCH --array=5-6
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_%a.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_%a.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Exact minimum number of monochromatic lines in 2-colourings of [3]^n, n = array index (5, 6).
# One single-threaded MaxSAT run per task; measured ~50 MB at n=5, so 1G each.
cd /home/hrmiller/CDCL-Experimentation
python3 -u experiments/hj_min_mono.py $SLURM_ARRAY_TASK_ID

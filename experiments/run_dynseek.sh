#!/bin/bash
#SBATCH --job-name=dynseek
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=9G
#SBATCH --time=01:30:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/dynseek.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/dynseek.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

cd /home/hrmiller/CDCL-Experimentation

# Dynamic "fail-first" variable choice (argmax S, ties by incidence order)
# with the avoid and saved value rules.
python3 -u experiments/hj_value_rule_experiment_v4.py \
  --instances 3:5,3:6,3:7,3:8,4:6,4:7 --orders dynseek --arms avoid,saved \
  --seeds 10 --workers 32 --timeout 1200 \
  --out experiments/results/dynseek.csv --resume

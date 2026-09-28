#!/bin/bash
#SBATCH --job-name=valrule
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/value_rule_overnight.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/value_rule_overnight.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

cd /home/hrmiller/CDCL-Experimentation
OUT=experiments/results/value_rule_overnight.csv
INST=4:5,4:6,3:5,3:6,3:7,3:8,4:7,4:8     # hardest (4:8) last

# Part 1: all five arms with the incidence order (the CaDiCaL arm ignores the
#         order, so it only needs to run once)
python3 -u experiments/hj_value_rule_experiment_v2.py --instances $INST \
  --orders incidence --arms cadical,constant,random,seek,avoid \
  --seeds 10 --workers 32 --timeout 1200 --out $OUT --resume

# Part 2: the four value rules under two other fixed variable orders
python3 -u experiments/hj_value_rule_experiment_v2.py --instances $INST \
  --orders index,random --arms constant,random,seek,avoid \
  --seeds 10 --workers 32 --timeout 1200 --out $OUT --resume

#!/bin/bash
#SBATCH --job-name=mmdec5
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=02:00:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n5_%a.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n5_%a.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Calibration on n=5 (known answer: at most 14 is UNSAT). Task 0: symmetry breaking on; task 1: off.
cd /home/hrmiller/CDCL-Experimentation
if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then SB=--sb; else SB=; fi
python3 -u experiments/hj_min_mono_dec.py --n 5 --bound 14 $SB --kissat $HOME/kissat/build/kissat

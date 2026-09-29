#!/bin/bash
#SBATCH --job-name=mmdec6
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500M
#SBATCH --time=48:00:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n6.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n6.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Is there a 2-colouring of [3]^6 with at most 134 monochromatic lines? (digit-sum colouring has 135)
# UNSAT => the digit-sum colouring is optimal at n=6. Kissat, symmetry breaking on.
# Measured: kissat ~164 MB after 45 s on this 441k-variable CNF; 1500M leaves room for growth.
cd /home/hrmiller/CDCL-Experimentation
python3 -u experiments/hj_min_mono_dec.py --n 6 --bound 134 --sb --kissat $HOME/kissat/build/kissat

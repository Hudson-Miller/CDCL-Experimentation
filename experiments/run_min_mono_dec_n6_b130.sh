#!/bin/bash
#SBATCH --job-name=mmdec6b130
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500M
#SBATCH --time=48:00:00
#SBATCH --output=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n6_b130.slurm.out
#SBATCH --error=/home/hrmiller/CDCL-Experimentation/experiments/results/min_mono_dec_n6_b130.slurm.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hrmiller@colgate.edu

# Is there a 2-colouring of [3]^6 with at most 130 monochromatic lines? (best symmetric colouring has 131)
# UNSAT => M(6) = 131 and the optimum is symmetric. Kissat, symmetry breaking on.
# Measured: kissat ~164 MB after 45 s on this 441k-variable CNF; 1500M leaves room for growth.
cd /home/hrmiller/CDCL-Experimentation
python3 -u experiments/hj_min_mono_dec.py --n 6 --bound 130 --sb --kissat $HOME/kissat/build/kissat

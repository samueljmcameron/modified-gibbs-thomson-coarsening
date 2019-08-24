#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out
#SBATCH --array=0-7

module restore standard_modules

python single-sigma.py $1 $SLURM_ARRAY_TASK_ID

#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A_%a.out
#SBATCH --array=0-32

module restore standard_modules

python single_sigma.py $1 $SLURM_ARRAY_TASK_ID
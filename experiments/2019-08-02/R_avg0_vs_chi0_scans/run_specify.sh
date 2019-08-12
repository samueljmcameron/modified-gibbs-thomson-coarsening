#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A.out

module restore standard_modules

python single_sigma_specify_chi_and_R0.py $1 $2

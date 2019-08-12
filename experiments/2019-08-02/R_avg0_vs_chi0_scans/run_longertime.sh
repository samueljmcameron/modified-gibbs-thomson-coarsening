#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-arutenbe
#SBATCH --output=slurmoutput/run_%A.out

module restore standard_modules

python single_sigma.py $1 $2

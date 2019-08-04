#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=def-arutenbe

module restore standard_modules

python organize_files.py $1

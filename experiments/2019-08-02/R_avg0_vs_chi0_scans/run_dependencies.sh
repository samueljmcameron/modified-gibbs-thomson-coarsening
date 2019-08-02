#!/bin/bash

mkdir -p slurmoutput

jout1=$(sbatch run_first.h $1)

jid1="${jout1//[!0-9]/}"

echo $jid1

jout2=$(sbatch --dependency=afterok:$jid1 run_second.h $1)

import os
from pathlib import Path
import numpy as np
import subprocess
import time

def find_chi_0(key,dir_name):

    with open(dir_name+f"run_{key}_0.out") as fp:

        linesto7 = [next(fp) for x in range(7)]


        line7 = linesto7[-1]

    return line7[6:]

def chi_index(key,dir_name):


    chi_0 = float(find_chi_0(key,dir_name))

    chi_0s = np.linspace(-10,25,num=36,endpoint=True)

    idx = np.abs(chi_0s-chi_0).argmin()

    return idx

def R_index(fname):

    a = fname[13:15]

    if a[-1] == '.':

        a = a[0]

    return int(a)

if __name__=="__main__":



    dir_name = str(Path.cwd()) + "/slurmoutput/"

    directory = os.fsencode(dir_name)


    bad_files = {}

    for file in os.listdir(directory):

        fname = os.fsdecode(file)

        key = fname[4:12]

        with open(dir_name+fname) as fp:

            line1=fp.readline()
            line2=fp.readline()

            if line2[:10] == "slurmstepd":

                subprocess.run(["sbatch","run_longertime.sh",str(chi_index(key,dir_name)),
                               str(R_index(fname))],check=True)
                bad_files[fname] = [chi_index(key,dir_name),R_index(fname)]
                time.sleep(0.5)

                

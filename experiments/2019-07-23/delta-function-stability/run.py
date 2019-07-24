import numpy as np
import sys

sys.path.append('../../scripts/')

from onerun import OneRun

def chi_0(R_avg0):

    R_eq = 10.0

    return (5/3*R_avg0-R_eq)*(R_avg0-R_eq)


if __name__=="__main__":

    R_avg0 = float(sys.argv[1])

    scan = {}

    scan['R_avg0'] = str(R_avg0)

    scan['chi_0'] = str(chi_0(R_avg0))

    run = OneRun(scan=scan)

    run.run_exe()

    run.mv_standard_files()

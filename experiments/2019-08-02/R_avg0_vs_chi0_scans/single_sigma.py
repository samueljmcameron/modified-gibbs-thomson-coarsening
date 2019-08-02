import numpy as np
import sys

sys.path.append('../../scripts/')

from onerun import OneRun


# look along stability where Rdot = 0 (in zero vol frac limit) and 
# see if making sigma bigger lets fixed point move.


if __name__=="__main__":

    R_avg0s = np.linspace(4,12,num=33,endpoint=True)

    chi_0s = np.linspace(-10,25,num=36,endpoint=True)

    chi_0 = chi_0s[int(sys.argv[1])]

    R_avg0 = R_avg0s[int(sys.argv[2])]

    scan = {}

    scan['R_avg0'] = str(R_avg0)

    scan['chi_0'] = str(chi_0)

    sigma = 0.01

    scan['sigma_R0'] = str(sigma)

    run = OneRun(scan=scan,executable = "../../../bin/3d-gaussian-ic-only-finaltime")

    # run calculation

    run.run_exe()

    # move output file from tmp folder to data folder

    fname="scanning"

    run.mv_file(fname)


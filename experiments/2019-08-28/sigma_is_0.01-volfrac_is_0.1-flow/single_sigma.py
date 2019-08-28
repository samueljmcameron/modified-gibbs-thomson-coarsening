import numpy as np
import sys

sys.path.append('../../scripts/')

from onerun import OneRun


# look along stability where Rdot = 0 (in zero vol frac limit) and 
# see if making sigma bigger lets fixed point move.


if __name__=="__main__":

    R_avg0 = float(sys.argv[1])

    index_chi = int(sys.argv[2])

    chi_0s = np.linspace(-10,25,num=8,endpoint=True)
    
    chi_0 = chi_0s[index_chi]

    scan = {}

    scan['R_avg0'] = str(R_avg0)

    scan['chi_0'] = str(chi_0)

    run = OneRun(scan=scan,executable="../../../bin/3d-gaussian-cubic-flow")

    # run calculation

    run.run_exe()

    run.mv_file("obs_vs_t")


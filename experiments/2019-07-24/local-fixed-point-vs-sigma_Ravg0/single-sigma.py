import numpy as np
import sys

sys.path.append('../../scripts/')

from onerun import OneRun

def chi_0(R_avg0):

    R_eq = 10.0

    return (5/3*R_avg0-R_eq)*(R_avg0-R_eq)


# look along stability where Rdot = 0 (in zero vol frac limit) and 
# see if making sigma bigger lets fixed point move.


if __name__=="__main__":

    R_avg0 = float(sys.argv[1])

    scan = {}

    scan['R_avg0'] = str(R_avg0)

    scan['chi_0'] = str(chi_0(R_avg0))

    sigmas = np.array([0,0.001,0.01,0.1,1.0],float)

    sigma = sigmas[int(sys.argv[2])]

    scan['sigma_Ravg0'] = str(sigma)

    run = OneRun(scan=scan)

    # run calculation

    run.run_exe()

    t_intervals = int(run.params['t_intervals'])

    # move chi vs t to data folder

    chi_name="chi_vs_t"

    run.mv_file(chi_name)

    # move initial radius distribution to data folder

    r0_name = "radius_0"

    run.mv_file(r0_name)

    # move last two saved radius distribtions to see if in steady state

    rsecondlast_name = f"radius_{t_intervals-1}"

    run.mv_file(rsecondlast_name)

    rfinal_name = f"radius_{t_intervals}"

    run.mv_file(rfinal_name)

    # move initial basis distribution to data folder

    b0_name = "basis_0"

    run.mv_file(b0_name)


    # move last two saved basis distributions

    bsecondlast_name = f"basis_{t_intervals-1}"

    run.mv_file(bsecondlast_name)

    bfinal_name = f"basis_{t_intervals}"

    run.mv_file(bfinal_name)

    # remove other files from the tmp_data folder.

    #run.delete_standard_files_from_tmp_path()

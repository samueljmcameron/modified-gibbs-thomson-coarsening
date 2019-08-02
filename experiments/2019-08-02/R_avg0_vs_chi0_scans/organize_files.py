import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

sys.path.append('../../scripts/')

from loaddata import LoadData,sort_R_avg0_vs_chi_0_data



if __name__=="__main__":

    sigma = 0.01

    savesuf = ["R_eq","volFrac_0","beta","chi_0"]

    R_avg0s = np.linspace(4,12,num=33,endpoint=True)

    chi_0s = np.linspace(-10,25,num=36,endpoint=True)

    chi_0 = chi_0s[int(sys.argv[1])]

    scan = {}

    scan['sigma_R0'] = str(sigma)

    scan['chi_0'] = str(chi)

    sort_R_avg0_vs_chi_0_data(R_avg0s,scan=scan,savesuf = ["R_eq","volFrac_0","beta","chi_0"])

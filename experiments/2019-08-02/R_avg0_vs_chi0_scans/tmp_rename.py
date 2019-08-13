import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

sys.path.append('../../scripts/')

from loaddata import tmp_rename_data



if __name__=="__main__":


    chi_0s = np.linspace(-10,25,num=36,endpoint=True)

    tmp_rename_data(chi_0s)

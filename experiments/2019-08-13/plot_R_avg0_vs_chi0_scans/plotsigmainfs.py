import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sys.path.append('../../scripts/')

from fig_settings import fig_settings
from loaddata import LoadData
from plotNinfs import quadratic, get_2d_array


if __name__=="__main__":


    fig_settings()

    colors = sns.color_palette("muted",8)

    width = 3.37
    height = 3.37

    fig, ax = plt.subplots()

    fig.set_size_inches(width,height)

    chi0s = np.linspace(-10,25,num=36,endpoint=True)

    R_avg0s = np.linspace(4,12,num=33,endpoint=True)

    RRs,ld = get_2d_array(3,chi0s,R_avg0s)

    cs = ax.contourf(R_avg0s,chi0s,RRs,levels=30)
    ax.plot(R_avg0s,quadratic(R_avg0s,float(ld.params['R_eq'])),'k-',lw=4)

    cbar = fig.colorbar(cs)

    ax.set_xlabel(r'$R(0)$')
    ax.set_ylabel(r'$\chi(0)$')

    cbar.ax.set_ylabel(r"$\sigma(t\to\infty)$")

    fig.subplots_adjust(left=0.18,bottom=0.15)
    fig.savefig(ld.file_savename("sigmainfty"))

    plt.show()

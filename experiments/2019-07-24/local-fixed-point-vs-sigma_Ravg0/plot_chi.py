import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from single_sigma import chi_0

sys.path.append('../../scripts/')

from fig_settings import fig_settings
from loaddata import LoadData

if __name__=="__main__":

    fig_settings()


    width = 16
    height = 10

    colors = sns.color_palette("muted",8)

    savesuf = ["R_eq","volFrac_0","beta"]


    R_avg0s = np.linspace(5,12,num=8,endpoint=True)

    sigmas = np.array([0,0.001,0.01,0.1,1.0],float)

    fig, axarr = plt.subplots(2,4)

    fig.set_size_inches(width,height)

    scan = {}

    lines = [":","-.","--","-",":"]

    for i,R_avg0 in enumerate(R_avg0s):

        scan['R_avg0'] = str(R_avg0)

        scan['chi_0'] = chi_0(R_avg0)

        for j,sigma in enumerate(sigmas):

            scan['sigma_Ravg0'] = str(sigma)
        
            ld = LoadData(scan = scan,savesuf=savesuf)

            ts = ld.data[:,0]

            chis = ld.data[:,1]

            axarr.flat[i].set_title(rf"$<\,R\,>(t)={R_avg0:.1f}$")
            axarr.flat[i].plot(ts,chis,'-',color=colors[j],linestyle=lines[j],
                               lw=4,
                               label=rf"$\sigma=\num{{{sigma:.0e}}}$")

        axarr.flat[i].set_xlabel(r"$t$")
        axarr.flat[i].set_ylabel(r"$\chi(t)$")

    axarr.flat[i].legend(frameon=False,handlelength=5)

    #    for ax in axarr.flat:
    #        ax.label_outer()

    fig.subplots_adjust(bottom = 0.08,top = 0.95,left=0.08,right=0.95)

    plt.show()

import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append('../../2019-07-24/local-fixed-point-vs-sigma_Ravg0/')

from single_sigma import chi_0

sys.path.append('../../scripts/')

from fig_settings import fig_settings
from loaddata import LoadData
from readparams import ReadParams

if __name__=="__main__":


    fig_settings()


    width = 16
    height = 10

    colors = sns.color_palette("muted",8)

    savesuf = ["R_eq","volFrac_0","beta"]

    data_path = "../../2019-07-24/local-fixed-point-vs-sigma_Ravg0"

    loadfilepath = data_path + "/data"

    datfile = data_path + "/data/input.dat"
    
    R_avg0s = np.linspace(5,12,num=8,endpoint=True)

    sigmas = np.array([0,0.001,0.01,0.1,1.0],float)

    fig, axarr = plt.subplots(2,4)

    fig.set_size_inches(width,height)

    scan = {}

    markers = ["o-","v-","^-","<-",">-"]

    for i,R_avg0 in enumerate(R_avg0s):

        scan['R_avg0'] = str(R_avg0)

        scan['chi_0'] = chi_0(R_avg0)



        for j,sigma in enumerate(sigmas):

            scan['sigma_Ravg0'] = str(sigma)
        
            rp = ReadParams(scan=scan,datfile=datfile)

            ts = rp.list_of_t_vals()

            sigma_Ravgts = np.empty([len(ts)],float)


            for i_t,t in enumerate(ts):

                ld = LoadData(name=f"radius_{i_t}",scan=scan,savesuf=savesuf,
                              loadfilepath=loadfilepath,datfile=datfile)
                
                Rs = ld.data[:,1]

                sigma_Ravgts[i_t] = Rs.std()

            axarr.flat[i].plot(ts,sigma_Ravgts,markers[j],color=colors[j],
                               label=rf"$\sigma=\num{{{sigma:.0e}}}$")            



    
            axarr.flat[i].set_title(rf"$<\,R\,>(t=0)={R_avg0:.1f}$")
            axarr.flat[i].set_xscale('log')
            axarr.flat[i].set_yscale('log')

        axarr.flat[i].set_xlabel(r"$t$")
        axarr.flat[i].set_ylabel(r"$\sigma(t)$")

    axarr.flat[i].legend(frameon=False,handlelength=5)

    #    for ax in axarr.flat:
    #        ax.label_outer()

    fig.subplots_adjust(bottom = 0.08,top = 0.95,left=0.08,right=0.95)

    fig.savefig(ld.file_savename("sigma_Ravg-loglog"))

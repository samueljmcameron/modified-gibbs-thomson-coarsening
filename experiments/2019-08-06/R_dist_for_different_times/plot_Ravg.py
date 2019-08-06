import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sys.path.append('../../2019-07-24/local-fixed-point-vs-sigma_R0/')

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

    data_path = "../../2019-07-24/local-fixed-point-vs-sigma_R0"

    loadfilepath = data_path + "/data"

    datfile = data_path + "/data/input.dat"

    R_avg0s = np.linspace(5,12,num=8,endpoint=True)

    sigmas = np.array([0,0.001,0.01,0.1,1.0],float)

    scan = {}


    markers = ["o-","v-","^-","<-",">-"]

    for i,R_avg0 in enumerate(R_avg0s):

        scan['R_avg0'] = str(R_avg0)

        scan['chi_0'] = str(chi_0(R_avg0))



        for j,sigma in enumerate(sigmas):

            scan['sigma_R0'] = str(sigma)
        
            rp = ReadParams(scan=scan,datfile=datfile)

            ts = rp.list_of_t_vals()
            
            print(len(ts))

            fig, axarr = plt.subplots(3,len(ts)//3)

            fig.set_size_inches(width,height)

            for i_t,t in enumerate(ts):

                ld = LoadData(name=f"radius_{i_t}",scan=scan,loadfilepath=loadfilepath,
                              datfile=datfile)
                
                Rs = ld.data[:,1]


                axarr.flat[i_t].hist(Rs,density=True,stacked=True)

                axarr.flat[i_t].set_xlabel(r"$<R>(t)$")
                axarr.flat[i_t].set_xlim(0,20)
                axarr.flat[i_t].set_ylim(0,1)
                axarr.flat[i_t].text(3,0.5,rf"$t=\num{{{t:.1e}}}$")

    
            fig.suptitle(rf"$<\,R\,>(t=0)={R_avg0:.1f}$"
                         +", "+rf"$\sigma(t=0)=\num{{{sigma:.0e}}}$")

            fig.subplots_adjust(bottom = 0.08,top = 0.95,left=0.08,right=0.95)

            fig.savefig(ld.file_savename("R_avg"))

            plt.close(fig=fig)

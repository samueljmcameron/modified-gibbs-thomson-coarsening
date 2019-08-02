import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

    savesuf = ["R_eq","volFrac_0","beta"]

    data_path = "../../2019-07-24/local-fixed-point-vs-sigma_R0"

    loadfilepath = data_path + "/data"

    datfile = data_path + "/data/input.dat"

    
    R_avg0s = np.linspace(5,12,num=8,endpoint=True)

    sigmas = np.array([0,0.001,0.01,0.1,1.0],float)

    #fig, axarr = plt.subplots(2,4)

    #fig.set_size_inches(width,height)

    fig = plt.figure()

    ax = fig.add_subplot(111,projection='3d')

    Rs = np.linspace(4,13,num=100,endpoint=True)

    ax.plot(Rs,1000*np.ones(len(Rs)),chi_0(Rs),'k-')


    scan = {}


    markers = ["o-","v-","^-","<-",">-"]

    for i,R_avg0 in enumerate(R_avg0s):

        scan['R_avg0'] = str(R_avg0)

        scan['chi_0'] = chi_0(R_avg0)



        for j,sigma in enumerate(sigmas):

            scan['sigma_R0'] = str(sigma)
        
            rp = ReadParams(scan=scan,datfile=datfile)

            ts = rp.list_of_t_vals()

            R_avgts = np.empty([len(ts)],float)
            N_ts = np.empty([len(ts)],float)

            chi_spaced = np.empty([len(ts)],float)
            
            ldchi = LoadData(scan=scan,savesuf=savesuf,loadfilepath=loadfilepath,
                             datfile=datfile)

            t_smalls = ldchi.data[:,0]
            chis = ldchi.data[:,1]


            for i_t,t in enumerate(ts):

                ld = LoadData(name=f"radius_{i_t}",scan=scan,savesuf=savesuf,
                              loadfilepath=loadfilepath,datfile=datfile)
                
                Rs = ld.data[:,1]

                R_avgts[i_t] = Rs.mean()

                N_ts[i_t] = len(Rs)
                
                chi_spaced[i_t] = chis[t_smalls<=t][-1]



            #axarr.flat[i].plot(chi_spaced,R_avgts,markers[j],color=colors[j],
            #                   label=rf"$\sigma=\num{{{sigma:.0e}}}$")            

            if i == 0:
                label = rf"$\sigma=\num{{{sigma:.0e}}}$"
            else:
                label = None


            ax.plot(R_avgts,N_ts,chi_spaced,lw=4,color=colors[j],
                    label=label)

    
            #axarr.flat[i].set_title(rf"$<\,R\,>(t=0)={R_avg0:.1f}$")            




    ax.set_xlabel(r"$<R>(t)$")
    ax.set_ylabel(r"$N(t)$")
    ax.set_zlabel(r"$\chi(t)$")

    ax.legend(frameon=False,handlelength=5)

    #    for ax in axarr.flat:
    #        ax.label_outer()
    
    ax.view_init(30,-45)


    fig.subplots_adjust(bottom = 0.08,top = 0.95,left=0.08,right=0.95)

    fig.savefig(ld.file_savename("R_avg-vs-N-vs-chi"))

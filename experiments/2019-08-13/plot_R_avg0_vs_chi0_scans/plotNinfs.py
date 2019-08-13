import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns

sys.path.append('../../scripts/')

from fig_settings import fig_settings
from loaddata import LoadData

def quadratic(R,Req):

    return (5/3*R-Req)*(R-Req)


def get_2d_array(index,chi0s,R_avg0s):


    data_path = "../../2019-08-02/R_avg0_vs_chi0_scans"

    loadfilepath = data_path + "/data_formatted"

    datfile = data_path + "/data_formatted/input.dat"

    scan = {}


    loadsuf=["sigma_R0","R_eq","volFrac_0","beta","chi_0"]
    savesuf=["sigma_R0","R_eq","volFrac_0","beta"]

    thearray = np.zeros([len(chi0s),len(R_avg0s)],float)

    for i,chi0 in enumerate(chi0s):

        scan['chi_0'] = str(chi0)

        ld = LoadData(name=f"scanning",scan=scan,loadfilepath=loadfilepath,
                      loadsuf=loadsuf,savesuf=savesuf,
                      datfile=datfile)

        row = ld.data[:,index]
        
        if index != 4 and index != -1:

            Nrow = ld.data[:,4].astype(int)

            row[np.where(Nrow <= 1)] = np.nan
        
        thearray[i,:] = row
        

    return thearray,ld

if __name__=="__main__":


    fig_settings()

    colors = sns.color_palette("muted",8)

    width = 3.37
    height = 3.37

    fig, ax = plt.subplots()

    fig.set_size_inches(width,height)

    chi0s = np.linspace(-10,25,num=36,endpoint=True)

    R_avg0s = np.linspace(4,12,num=33,endpoint=True)

    NNs,ld = get_2d_array(4,chi0s,R_avg0s)

    cs = ax.contourf(R_avg0s,chi0s,NNs,levels=30)
    ax.plot(R_avg0s,quadratic(R_avg0s,float(ld.params['R_eq'])),'k-',lw=4)

    cbar = fig.colorbar(cs)

    ax.set_xlabel(r'$R(0)$')
    ax.set_ylabel(r'$\chi(0)$')

    cbar.ax.set_ylabel(r"$N(t\to\infty)$")

    fig.subplots_adjust(left=0.18,bottom=0.15)
    fig.savefig(ld.file_savename("Ninfty"))

    plt.show()

import numpy as np
import subprocess
import sys
import os
from readparams import ReadParams

class LoadData(ReadParams):
    # Perform a single run of calculation for the radius distribution.
    # You can choose which parameters to vary by specifying a dictionary
    # of values "scan" when creating the object.

    # attributes:
    #  datfile - the data file that you read parameters from
    #  scan - a dictionary where you can pre-specify certain parameters (vs
    #         reading them in through the datfile).
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,scan={},
                 datfile="data/input.dat",
                 loadsuf=["R_avg0","sigma_R0","R_eq","volFrac_0","beta",
                          "chi_0"],
                 savesuf=["R_avg0","sigma_R0","R_eq","volFrac_0","beta",
                          "chi_0"],
                 name="chi_vs_t",loadfilepath="data",savefilepath="results"):

        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.name = name
        self.loadfilepath = loadfilepath
        self.savefilepath = savefilepath
        if os.path.isfile(self.file_name()):
            self.data = np.loadtxt(self.file_name())
            self.overlap=False
        elif os.path.isfile(self.file_name(overlap=True)):
            self.data = np.loadtxt(self.file_name(overlap=True))
            print(f"Found file {self.file_name()}, but drops are overlapping!")
            self.overlap=True
        else:
            print(f"Could not find file {self.file_name()}")
            self.overlap=False

        return
        

    def file_name(self,overlap=False):
    
        suffix = self.write_suffix()

        if overlap:
            suffix = suffix + "-overlap"


        fname =f"{self.loadfilepath}/{self.name}_{suffix}.txt"

        return fname


    def file_savename(self,varname,file_format="pdf"):

        suffix = self.write_suffix(suffix_type="save")

        return f"{self.savefilepath}/{varname}_{suffix}.{file_format}"

    def remove_file(self,overlap=False):

        os.remove(self.file_name(overlap=overlap))

        return

    def read_header(self,overlap=False):

        with open(self.file_name(overlap=overlap)) as f:
            first_line = f.readline()

        return first_line





def sort_R_avg0_vs_chi_0_data(R_avg0s,scan={},savesuf=["R_eq","volFrac_0","beta","chi_0"]):

    # create new array to store all of the (maybe) converged data

    new_data = np.empty([len(R_avg0s),7],float)

    for i,R_avg in enumerate(R_avg0s):

        scan['R_avg0'] = str(R_avg)

        ld = LoadData(name="scanning",scan=scan,savefilepath="data",savesuf=savesuf)

        if i == 0:

            header = ld.read_header()

        if ld.data.ndim == 1:

            print("potentially didn't converge, not enough data to tell!")

            time = ld.data

            convergence_flag = 2

        else:

            almost_time = ld.data[-2,:]

            time = ld.data[-1,:]

            if np.any(np.abs(almost_time-time)>1e-15):

                convergence_flag = 1

            else:

                convergence_flag = 0

                ld.remove_file()


        new_data[i,:] = np.concatenate((time,[convergence_flag]))

    np.savetxt(ld.file_savename("scanning",file_format="txt"),new_data,
               fmt="%e\t%e\t%e\t%e\t%e\t%e\t%d",header=header)

    return

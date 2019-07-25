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
                 loadsuf=["R_avg0","sigma_Ravg0","R_eq","volFrac_0","beta",
                          "chi_0"],
                 savesuf=["R_avg0","sigma_Ravg0","R_eq","volFrac_0","beta",
                          "chi_0"],
                 name="chi_vs_t",loadfilepath="data",savefilepath="results"):

        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.name = name
        self.loadfilepath = loadfilepath
        self.savefilepath = savefilepath
        if os.path.isfile(self.file_name()):
            self.data = np.loadtxt(self.file_name())
        else:
            print(f"Could not find file {self.file_name()}")

        return
        

    def file_name(self):
    
        suffix = self.write_suffix()

        fname =f"{self.loadfilepath}/{self.name}_{suffix}.txt"

        return fname


    def file_savename(self,varname,file_format="pdf"):

        suffix = self.write_suffix(suffix_type="save")

        return f"{self.savefilepath}/{varname}_{suffix}.{file_format}"

import numpy as np


class ReadParams(object):

    def __init__(self,datfile="data/input.dat",scan={},
                 loadsuf=["R_avg0","sigma_Ravg0","R_eq","volFrac_0","beta",
                          "chi_0","N","dt_max","t_final","t_intervals"],
                 savesuf=["R_avg0","sigma_Ravg0","R_eq","volFrac_0","beta",
                          "chi_0","N","dt_max","t_final","t_intervals"]):

        self.datfile = datfile
        self.scan = scan
        self.loadsuf = loadsuf
        self.savesuf = savesuf
        self.params = self.read_params()
        
        return

    def set_param(self,key,val):
        # if the key passed as an argument here is already a key in the
        # self.scan dictionary, return the value in self.scan[key].
        # Otherwise, return val.
        
        for scan_key,scan_val in self.scan.items():

            if (key == scan_key):
                return scan_val

        return val
        

    def read_params(self):
        # read input parameters from datfile, unless the input parameter
        # is specified in self.scan dictionary, in which case, read it
        # from there. return a dictionary with parameter labels and names.
        # NOTE since this is running python 3.7 (or later), dicts are ordered.
        
        params = {}
        
        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()

                params[key]=self.set_param(key,val)

        return params

    def write_suffix(self,suffix_type="load"):
        # write the ending "suffix" of the file name (all of the trailing parameter values
        # in the filename. e.g. if fname = "energy_3.00000e+00_5.32000e-01.txt", then the
        # suffix="3.00000e+00_5.32000e-01", WITHOUT the ".txt".

        if (suffix_type=="save"):
            cpylist = self.savesuf
        else:
            cpylist = self.loadsuf

        plist = [float(self.params[s]) for s in cpylist]

        suffix = '_'.join([f'{float(s):.4e}' for s in plist])
        
        return suffix

    def list_of_t_vals(self,second_t_val = 1.0):

        t0 = second_t_val

        t_len = int(self.params['t_intervals'])+1

        t_final = float(self.params['t_final'])

        ts = np.empty([t_len],float)

        ts[0] = 0.0

        exp10_t0 = np.log10(second_t_val)
        
        exp10_tf = np.log10(t_final)
        
        ts[1:] = np.logspace(0,np.log10(t_final),num=t_len-1,endpoint=True)

        return ts

import extract_spectral_data
import parameter_ranges as pr
import numpy as np
import pandas as pd
from scipy import interpolate
import time
import sys

print("imports complete")

class RunData:
    def __init__(self,n,i,t) :
        """run_setup(density,intensity,temperature)

        Calculate the directories, and extract radii & column depths for a given run
        """
        self.n=n
        self.i=i
        self.t=t

        self.data_dir = pr.format_dir(self.n,self.i,self.t)
        self.spec_file_base = self.data_dir + pr.format_file(self.n,self.i,self.t) + ".znu"

        # extract depths and distances from one file (should all be identical!)
        file = f"{self.spec_file_base}0"
        self.radii,self.depths = np.loadtxt(file, usecols=[0, 1], skiprows=1,unpack=True)


def extract_heat_cool(rd):
    cool_file = rd.data_dir + pr.format_file(rd) + ".cool"
    return pd.read_csv(cool_file,sep='\t',usecols=[2,3],names=['Htot','Cool'],skiprows=[0],comment='#')

def extract_dust(rd):
    dg_filename = rd.data_dir + pr.format_file(rd) + ".ratio"
    dg_matrix=pd.read_csv(dg_filename,sep='\t',usecols=range(1,22),comment='#',header=None).values

    Tdust_filename = rd.data_dir + pr.format_file(rd) + ".gtemp_full"
    Tdust_matrix = pd.read_csv(Tdust_filename,sep='\t',usecols=range(1,21),comment='#',header=None).values

    # if mass not taken into account, weight mean area by <a>**3/<a>**2 i.e. mass
    # we use dg_matrix for mass, so we don't need to do this
    dust_area_weights = np.tile(pr.d_rads,2)
    dust_area_weights = dust_area_weights/np.sum(dust_area_weights)

    # weight by abundance and emission
    Tdust = np.sum(Tdust_matrix**4*dust_area_weights*dg_matrix[:,:-1],axis=1)**.25
    dg = dg_matrix[:,-1]
    return pd.DataFrame(np.array([Tdust,dg]).T,columns=["Tdust","dg"])

def extra_one_run(n,i,t):
    rd = RunData(n,i,t)

    spec = extract_spectral_data.SpectrumInterpolator(rd)
    spec.calculate_rad_pressure() #calculates all other spectral terms too

    df = spec.gen_dataframe()

    df = df.join(extract_heat_cool(rd))
    df = df.join(extract_dust(rd))

    df_out = pd.DataFrame()
    for key in df:
        df_out[key] = interpolate.interp1d(df["tau"],df[key])(pr.table_tau)
    return df_out

if __name__ == "__main__" :
    t0 = time.time()
    df_out = extra_one_run(1.,10.2,1.)
    print(time.time()-t0)
    df_out.to_csv("data/test.txt",sep=" ",index=False)


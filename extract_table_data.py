import extract_spectral_data
from parameters import parameter_ranges as pr
import numpy as np
import pandas as pd
from scipy import interpolate
import time
import multiprocessing

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

def extract_one_run(n,i,t):
    print(n,i,t)
    rd = RunData(n,i,t)

    spec = extract_spectral_data.SpectrumInterpolator(rd)
    spec.calculate_rad_pressure() #calculates all other spectral terms too

    df = spec.gen_dataframe()

    df = df.join(extract_heat_cool(rd))
    df = df.join(extract_dust(rd))

    df_out = pd.DataFrame()
    for key in df:
        if np.any(df[key]<0):
            print("<0 in ",n,i,t,key)
        df_out[key] = interpolate.interp1d(df["tau"],np.log10(df[key]),fill_value="extrapolate")(pr.table_tau)
    df_out["tau"] = pr.table_tau
    return n,i,t,df_out

def full_run(alln,alli,allt,key_order=None,Ncores=4):
    if key_order is None:
        key_order = pr.key_order

    n0, i0, t0 = np.meshgrid(alln,alli,allt)

    points = np.array([n0.ravel(), i0.ravel(), t0.ravel()]).T

    t0 = time.time()

    pool = multiprocessing.Pool(Ncores)

    df_list = pool.starmap(extract_one_run, points)

    full_df = pd.DataFrame()
    print(time.time() - t0)

    nit_table = []

    for n,i,t,df in df_list:
        df["n"] = np.full(len(df),n)
        df["i"] = np.full(len(df),i)
        df["t"] = np.full(len(df),t)
        full_df = full_df.append(df)


    with open("data/table_key.txt",'w') as keyf:
        keyf.write(f"{len(pr.table_tau)} {len(alln)} {len(alli)} {len(allt)}\n")
        for t in [pr.table_tau,alln,alli,allt]:
            for i in t:
                keyf.write(f"{i}\n")

    full_df.sort_values(by=["n","i","t"],inplace=True)
    full_df = full_df[key_order]

    print(time.time() - t0)
    full_df.to_csv("data/full_test.txt",sep=" ",index=False)
    print(time.time() - t0)

def test_run():
    t0 = time.time()
    df_out = extract_one_run(pr.logn[1], pr.logi[0], pr.logt[0])
    print(time.time() - t0)
    df_out.to_csv("data/test.txt",sep=" ",index=False)

if __name__ == "__main__" :
    # test_run()
    full_run(pr.logn[0 :2], pr.logi[0 :2], pr.logt[0 :2],Ncores=12) # truncated for test

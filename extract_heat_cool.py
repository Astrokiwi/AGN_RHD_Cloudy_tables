import pandas as pd
import numpy as np
import parameter_ranges as pr

def extract(rd):
    cool_file = rd.data_dir + pr.format_file(rd) + ".cool"
    return pd.read_csv(cool_file,sep='\t',usecols=[2,3],names=['Htot','Cool'],skiprows=[0],comment='#')

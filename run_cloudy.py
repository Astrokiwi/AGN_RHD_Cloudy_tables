#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 09:53:53 2017

@author: mv3e16
"""
import numpy as np
import os
import psutil
import parameter_ranges as pr
import multiprocessing
import subprocess
import sys

#%%
def grep():
    
    PROCNAME = "cloudy.exe"

    c=0
    for proc in psutil.process_iter():
        if proc.name() == PROCNAME:
            c+=1
    return float(c)
#%%
def runcloudy(point):
    
    n,i,t = point

    print(n,i,t)
    
    os.chdir(f"{pr.directory}/n{n:.0f}/In{i:.1f}/Te{t:.2f}")
    p = subprocess.run(["cloudy",f"n{n:.0f}_In{i:.1f}_Te{t:.2f}.in"],text=True,stderr=subprocess.STDOUT)
    # p = subprocess.run(["ls"], capture_output=True,shell=True,text=True) # shell=True for passing a string
    # p = subprocess.run(["pwd"], capture_output=True,shell=True,text=True) # shell=True for passing a string
    # p = subprocess.run(["echo hi; sleep 2; echo ho"], capture_output=True,shell=True,text=True) # shell=True for passing a string

    print(p.stdout)

    # os.system("cloudy *in &")


#%%
if __name__ == "__main__" :
    if len(sys.argv)>1:
        Ncores = int(sys.argv[1])
    else:
        Ncores = 4

    n0, i0, t0 = np.meshgrid(pr.logn[0:2], pr.logi[0:2], pr.logt[0:2])
    # n0, i0, t0 = np.meshgrid(pr.logn, pr.logi, pr.logt)

    points = np.array([n0.ravel(), i0.ravel(), t0.ravel()]).T

    points=points[:]
#%%

    pool = multiprocessing.Pool(Ncores)
    
    PoolReturn = pool.map(runcloudy,points)

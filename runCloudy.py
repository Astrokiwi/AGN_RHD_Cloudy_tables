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

Ncores = 4

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
    p = subprocess.run(["cloudy","*in"])
    # p = subprocess.run(["echo hi; sleep 2; echo ho"], capture_output=True,shell=True,text=True) # shell=True for passing a string

    # print(p.stdout)

    # os.system("cloudy *in &")


#%%
if __name__ == "__main__" :
    n0, i0, t0 = np.meshgrid(pr.logn, pr.logi, pr.logt)

    points = np.array([n0.ravel(), i0.ravel(), t0.ravel()]).T

    points=points[:]
#%%

    pool = multiprocessing.Pool(Ncores)
    
    PoolReturn = pool.map(runcloudy,points)

#%%
# c=0
# n=2
# t=logt[16]
# for i in logi:
#        # for t in logt[10:15]:
#             os.chdir("{0}n{1}/In{2:.1f}/Te{3:.2f}".format(directory,n,i,t))
#             c=c+1
#             os.system("cloudy *in &")

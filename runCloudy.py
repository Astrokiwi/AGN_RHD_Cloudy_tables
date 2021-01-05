#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 09:53:53 2017

@author: mv3e16
"""
import numpy as np
import os
import psutil

#%%
directory= '/home/mv3e16/filestore/cloudy_data/'

I0=5.6*10**7
logi = []
logi.append(np.log10(I0*10**1.5))
logi.append(np.log10(I0*10))
logi.append(np.log10(I0*10**0.5))

for i in range(14):
    a = I0 * 10**(-i*0.5)
    logi.append(np.log10(a))
    
logt=[]
T0=10

for i in range(17):
    T=T0*10**(i*0.25)
    logt.append(np.log10(T))        

for i in range(3):
    T = 10**(6+i)
    logt.append(np.log10(T))
    
    
densities=np.arange(0,8)


n0,i0,t0= np.meshgrid(densities,logi,logt)

points = np.array([n0.ravel(),i0.ravel(),t0.ravel()]).T
#%%
def grep():
    
    PROCNAME = "cloudy.exe"

    c=0
    for proc in psutil.process_iter():
        if proc.name() == PROCNAME:
            c=c+1
    return(float(c))
#%%
def runcloudy(point):
    
    n=point[0]
    i=point[1]
    t=point[2]
    
    os.chdir("{0}n{1:.0f}/In{2:.1f}/Te{3:.2f}".format(directory,n,i,t))
    os.system("cloudy *in &")
    
    return()


#%%
points= points[:]
#%%
'''
if __name__ == "__main__":
    
    import multiprocessing 
    
    Ncores = 3
    pool = multiprocessing.Pool(Ncores)
    
    PoolReturn = pool.map(runcloudy,points)
    '''
#%%
c=0
n=2
t=logt[16]
for i in logi:
       # for t in logt[10:15]:
            os.chdir("{0}n{1}/In{2:.1f}/Te{3:.2f}".format(directory,n,i,t))
            c=c+1
            os.system("cloudy *in &")
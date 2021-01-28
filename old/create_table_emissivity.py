#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
import matplotlib as mpl
mpl.use('Agg')
import pylab as P
import os 
import glob

directory = "/Volumes/marta_filestore-1/cloudy_data/"  
#%%

ryd2microns=13.6*1.24

#%%
def makedfnew(df,n,i,t):
    
    return(df.query('tgas=={0} and i0=={1} and n0=={2}'.format(t,i,n)))

#%%    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

    return idx

#%%
I0=5.6*10**7

logi = []
rsub = []
logi.append(np.log10(I0*10**1.5))
logi.append(np.log10(I0*10))
logi.append(np.log10(I0*10**0.5))

for i in range(14):
    a = I0 * 10**(-i*0.5)
    r = '{0:.2f}'.format(10**(i*0.5/2))
    logi.append(np.log10(a))
    rsub.append(r)
    
logt=[]

T0=10

for i in range(17):
    T=T0*10**(i*0.25)
    logt.append(np.log10(T))
    
for i in range(3):
    T = 10**(6+i)
    logt.append(np.log10(T))

#%%
def f(x,x_points,y_points):
    
     c=interpolate.interp1d(x_points, y_points, kind='nearest',fill_value='extrapolate')
    
     return c(x)
#%%
#sample densities intensity temperatures
n=0
i=logi[0]
t=8

#%%

#We want the emissivity at 12,18,850 microns. We need to find the corresponding Cloudy values.
# I found the frequency bins array by Cloudy and retrieved the closest values to the frequencies we want

nu_micron  = np.array([12,18,850])
nu_ryd = nu_micron / ryd2microns
#%%
a=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,0),sep='\t',usecols=[1,2,3,4])
frequency_rydberg = float(a.columns[-1][12:-4]) # this returns the value of the frequency in rydberg

#%%

directory_files = '{0}n{1}/In{2:.1f}/Te{3:.2f}/'.format(directory,n,i,t)
os.chdir(directory_files)
files = glob.glob('*znu*')

#%%len)

#all_nuryd = []


for j in files: 
    a=pd.read_csv(j,sep='\t',usecols=[1,2,3,4])
    frequency_rydberg = float(a.columns[-1][12:-4])

    all_nuryd.append(frequency_rydberg)
    
#I saved it as a text file at some point
#%%

b=np.loadtxt(directory+"all_nuryd.txt")
#%%

#we want 8 microns too (instead of 18. there was a misunderstanding)
file4 = files[find_nearest(b,8/ryd2microns)]

#%%
# FINDING the corresponding frequency files 

file1 = files[find_nearest(b,nu_ryd[0])] #@ 7.120164e-01 Ryd *.znu3616
file2 = files[find_nearest(b,nu_ryd[1])] #@ 1.067472e+00 Ryd *.znu3697
file3 = files[find_nearest(b,nu_ryd[2])]  # 5.038856e+01 Ryd *.znu4468
file4 = files[find_nearest(b,8/ryd2microns)] # *.znu3535

    
#%%
colden=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,2,7.7,1),sep='\t',usecols=[1],names=['colden'],comment='#')
colden=colden.drop_duplicates('colden').reset_index(drop=True)

l=len(colden)

#%%
def emissivity(n,i,t):

    colden= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[0,1],names=['depth','colden'],comment='#')
    
    nujnu12 = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,3616),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values
    #nujnu18 =  pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,3697),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values
    nujnu8 =  pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,3535),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values
    nujnu850 = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,4468),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values

    table = np.c_[colden.colden,nujnu12,nujnu8,nujnu850]

    df=pd.DataFrame(table,columns=['colden','nujnu12','nujnu8','nujnu850'])
    return(df)
#%%
colden=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,2,7.7,1),sep='\t',usecols=[1],names=['colden'],comment='#')
colden=colden.drop_duplicates('colden').reset_index(drop=True)                 
#%%
def df_interpolated(n,i,t):
    
    l=len(colden)
    
    df=emissivity(n,i,t)
    
    tot=[]

    for col in df.columns[1:]:
        tot.append(f(colden,df.colden,df[col]))
    
    tot=np.c_[tot].T[0]
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot],columns=['n0','i0','tgas','colden','nujnu12','nujnu8','nujnu850'])
    
    return(d)
#%%
def df_interpolated_all(n):
    a=[]
    for i in logi:
        for t in logt:
            a.append(df_interpolated(n,i,t))
            
    d=pd.DataFrame(np.vstack(a),columns=['n0','i0','tgas','colden','nujnu12','nujnu8','nujnu850'])
    return(d)
    
#%%
df0=df_interpolated_all(0)
print('n{} done'.format(0))

df1=df_interpolated_all(1)
print('n{} done'.format(1))

df2 = df_interpolated_all(2)
print('n{} done'.format(2))

df3=df_interpolated_all(3)
print('n{} done'.format(3))

df4=df_interpolated_all(4)
print('n{} done'.format(4))

df5=df_interpolated_all(5)
print('n{} done'.format(5))

df6=df_interpolated_all(6)
print('n{} done'.format(6))

df7=df_interpolated_all(7)
print('n{} done'.format(7))

#%%
a= [df0,df1,df2,df3,df4,df5,df6,df7]
d=pd.concat(a)
file= open('{0}emissivity_170620.txt'.format('/Users/marta/Desktop/'),"w")
file.write(d.to_string())
file.close()#%
#%%

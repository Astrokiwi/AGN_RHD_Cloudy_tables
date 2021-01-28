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

directory = "/home/mv3e16/filestore/cloudy_data/"  

#%%
def makedfnew(df,n,i,t):
    
    return(df.query('tgas=={0} and i0=={1} and n0=={2}'.format(t,i,n)))

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


sizes_in_cloudy = np.array([6.080e-03,8.991e-03,1.330e-02,1.966e-02,2.907e-02,4.299e-02,6.358e-02,9.402e-02,1.390e-01,2.056e-01])
delta = np.log10(0.25)-np.log10(0.005)
bins=delta/10
all_sizes = np.asarray([10**(np.log10(0.005)+bins*i*0.5) for i in range(0,21)])

delta_a = np.diff(all_sizes[range(0,21,2)])

mid_sizes =np.setdiff1d(all_sizes,all_sizes[(range(0,21,2))]) #these are the ten sizes of grain in cloudy

sizes_weights = np.power(mid_sizes,-1.5)*delta_a

#mrn_norm = 24.284271247461902 #int a^(-3.5) * a^(2) da from 0.005 to 0.25

mrn_norm=np.sum(sizes_weights) #24.400502799395426

sizes_weights2=np.concatenate([sizes_weights,sizes_weights])
mrn_norm2=2*mrn_norm
#%%
def cooling(n,i,t):

    Tgrain = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.gtemp_full'.format(directory,n,i,t),sep='\t',usecols=range(1,21),comment='#',header=None)
    Tgrain_depth = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.gtemp_full'.format(directory,n,i,t),sep='\t',usecols=[0],comment='#',names=['depth'])

    colden= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[0,1],names=['depth','colden'],comment='#')
    dustToGasRatio = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.ratio'.format(directory,n,i,t),sep='\t',usecols=range(1,21),comment='#',header=None)
    cool = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.cool'.format(directory,n,i,t),sep='\t',usecols=[0,1,2,3],names=['depth','Temp','Htot','Cool'],skiprows=[0],comment='#')

    kabs_I=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.kabs2'.format(directory,n,i,t),names=['kabs'],delimiter=r"\s+").as_matrix().ravel()
    kscat_I=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.kscat2'.format(directory,n,i,t),names=['kscat'],delimiter=r"\s+").as_matrix().ravel()
    prad=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.prad2'.format(directory,n,i,t),names=['kscat'],delimiter=r"\s+").as_matrix().ravel()
    tau_david=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.tau_david2'.format(directory,n,i,t),names=['kscat'],delimiter=r"\s+").as_matrix().ravel()

    
    dgtot = dustToGasRatio.cumsum(axis=1).iloc[:,-1]
    dgtot=dgtot.as_matrix()
    
    T4 = Tgrain**4
    tgtot = T4*sizes_weights2
    tgtot = tgtot.cumsum(axis=1).iloc[:,-1]
    tgtot = tgtot/mrn_norm2
    tgtot = tgtot**(0.25)
    tgtot = tgtot.as_matrix()

    #problema per n=0,In=4.7,Te=1 : tutte Nzone diverse
    if cool.shape[0]!= colden.shape[0]:
        
        Cool=f(colden.depth,cool.depth,cool.Cool)
        Heat=f(colden.depth,cool.depth,cool.Htot)
        tgtot=f(colden.depth,Tgrain_depth,tgtot)            
        dgtot=f(colden.depth,df4_depth,dgtot)
        table = np.c_[colden.colden,tgtot,Heat,Cool,prad,dgtot,kabs_I,kscat_I,tau_david]
        print('size not matching {0}_{1:2f}_{2}'.format(n,i,t))
        table=np.c_[colden.colden,tgtot,Heat,Cool,prad,dgtot,kabs_I,kscat_I,tau_david]
     
    else:
        table = np.c_[colden.colden,tgtot,cool.Htot,cool.Cool,prad,dgtot,kabs_I,kscat_I,tau_david]
        #table=np.c_[colden.colden,cool.Htot,cool.Cool,prad,kabs_I,kscat_I,tau_david]

    df=pd.DataFrame(table,columns=['colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
    #df=pd.DataFrame(table,columns=['colden','heat','cool','prad','kabs','kscat','tau'])
    df.loc[df['tau']<0,'tau']= f(df.query('tau<0').colden,df.query('tau>0').colden,df.query('tau>0').tau)# the first tau values were negative
    return(df)
#%%
#sample colden
colden=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format("/home/mv3e16/filestore/cloudy_data/",2,7.7,1),sep='\t',usecols=[1],names=['colden'],comment='#')
colden=colden.drop_duplicates('colden').reset_index(drop=True)                 
#%%
def df_interpolated(n,i,t):
    
    l=len(colden)
    
    df=cooling(n,i,t)
    
    tot=[]

    for col in df.columns[1:]:
        tot.append(f(colden,df.colden,df[col]))
    
    tot=np.c_[tot].T[0]
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot],columns=['n0','i0','tgas','colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
   

    return(d)
#%%
def df_interpolated_all(n):
    a=[]
    for i in logi:
        for t in logt:
            a.append(df_interpolated(n,i,t))
            
    d=pd.DataFrame(np.vstack(a),columns=['n0','i0','tgas','colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
    return(d)

#%%
l=len(colden)    


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

file= open('{0}100818.txt'.format(directory),"w")
file.write(d.to_string())
file.close()#%
#%%

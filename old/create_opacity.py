#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 09:53:53 2017

@author: mv3e16
"""
import numpy as np
import pandas as pd
import psutil
from scipy import interpolate
import matplotlib.pyplot as plt
import os
#%%
directory='/home/mv3e16/filestore/cloudy_data/'

#%%
light_speed_micron = 299792458000000 #micron/s
light_speed_cm = 2.997925*10**10 #cm/s
m_proton = 1.672621898*10**(-24)  # proton mass g
thomson = 6.652458715810**(-29)  # thomson cross section  m^2
sigma_stefan_boltz = 5.670519*10**(-5) # erg cm2 K-4 s-1
h_planck = 6.62607004*10**(-27) #erg s
kappa_boltz = 1.38064852*10**(-16) # erg K
pi=np.pi
ryd2nu=3.28*15
ryd2microns=13.6*1.24
BB_const = (2 * h_planck *(light_speed_micron)**2)*10**(+8) # erg s-1 cm-2 * E+8 (CONVERSION FACTOR)

#%% interpolation function
def f(x,x_points,y_points):
    
     c=interpolate.interp1d(x_points, y_points, kind='nearest',fill_value='extrapolate')
    
     return c(x)
#%%
def blackbody_lambda(A,T): #remember to give 10**T !!!!!!! 
    
    a=h_planck*light_speed_micron/(A*kappa_boltz*T)
    
    BB = BB_const/(A**5*(np.exp(a)-1))
    
    return(BB)
#%%    DEFINING INTESITY AND TEMPERATURE FIXED PARAMETERS
I0=5.6*10**7

logi = []
rsub = []

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
    
logt=logt[-4:]
    
#%%   
def grep():
    
    PROCNAME = "python"

    c=0
    for proc in psutil.process_iter():
        if proc.name() == PROCNAME:
            c=c+1
    return(c)


#%%   Taking the frequency array from cloudy c17
n=0
i=logi[0]
t=8
densities =np.arange(0,8)
#anu=pd.read_csv('{0}mu_values.csv'.format(directory),sep='\t',comment='#',names=['nu'],skiprows=np.arange(3,5277).tolist())
anu=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.con'.format('/home/mv3e16/filestore/cloudy_data/',n,i,t),sep='\t',usecols=[0],names=['nu'],comment='#')
deltanu=anu.diff()
deltanu=deltanu.fillna(deltanu.loc[1])
#%%
def interp_I0(n,i,t):
    
    df=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.con'.format(directory,n,i,t),sep='\t',usecols=[0,1],names=['nu','I0'],comment='#')
    If=f(anu,df.nu,df.I0)
    
    return(If.ravel())
#%%
def kappa_cloudy(n,i,t):

    
    kabs=[]
    kscat=[]
    table1=[]
    table2=[]
    
    
    for anui in range(len(anu)):
            df=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,anui),sep='\t',usecols=[3,4],names=['kabs','kscat'],skiprows=[0],comment='#')
                
                kabs=pd.DataFrame(df.kabs.as_matrix(),columns=['nu{}'.format(anui)])
                kscat=pd.DataFrame(df.kscat.as_matrix(),columns=['nu{}'.format(anui)])
                
                table1.append(kabs)
                table2.append(kscat)
                
    table1=pd.concat(table1,axis=1) 
    table2=pd.concat(table2,axis=1) 
    
    return(table1,table2) # !!!!!!!!!!!!!!! dimensions still to correct !!!!!!!!!!!!!!!

#%%
def tau_matrix(n,i,t): # returning a matrix of size (len(depth) rows X 5277 columns)
    
    #cloudy radius
    dr = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.dr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['depth'],comment='#')
    #cleansing
    dr=dr.drop_duplicates()
    dr.loc[-1]=dr.iloc[0][0]
    dr.index=dr.index+1
    dr=dr.sort_index().as_matrix()
    
    deltar=np.diff(dr[:,0])
    
    #checking if the cloudy simulation has converged
   df=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,0),sep='\t',usecols=[1,3,4],names=['depth','kabs','kscat'],skiprows=[0],comment='#')
    
    if df.shape[0]!=deltar.shape[0]:
        df=df.sort_values('depth')
        df=df.reset_index(drop=True)
        print('check the file n{0}_In{1:.1f}_Te{2:.2f}.dr'.format(n,i,t))
        deltar=df.depth.diff()
        deltar=deltar.fillna(deltar.loc[1]).as_matrix()
    
    a=[]

    for anui in range(len(anu)):  #retrieving the abs/scat opacities for each frequency
        df=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,anui),sep='\t',usecols=[1,3,4],names=['depth','kabs','kscat'],skiprows=[0],comment='#')
        
        df=df.sort_values('depth')
        df=df.reset_index(drop=True)
        
        k_abs=df.kabs.as_matrix()
        k_scat=df.kscat.as_matrix()
    
        #summing Ktot for the corresponding dr 
        tau = np.multiply((k_abs+k_scat),deltar)
    
        #put the result in an array with lenght = Nzones for each nu
        dftau= pd.DataFrame(tau,columns=['nu{}'.format(anui)])
        a.append(dftau)
    
    table=pd.concat(a,axis=1) # returning a matrix of size (len(depth) rows X 5277 columns (Ktot*DELTAR) .
    optical_depth=table.cumsum(axis=0) #integral performed. Result is the optical depth for each zone and each nu (Nzones rows X 5277 columns) 

    return(optical_depth)
#%%    
def continuum(n,i,t):
    
    incident=interp_I0(n,i,t)
    tau=tau_matrix(n,i,t)
    continuum_matrix=incident*np.exp(-tau) #return a matrix with (Nzone rows X 5277 columns)
    optical_depth=(-1)*np.multiply(tau,deltanu.as_matrix().ravel()).cumsum(axis=1).iloc[:,-1]  
    return(continuum_matrix,optical_depth)
#%%
 #returning an array Nzone X 5277 for Kabs and Kscat
 
def kappa_I(n,i,t):
    
    con=continuum(n,i,t)
    I_full=con[0]
    optical_depth=con[1]

    kcloudy=kappa_cloudy(n,i,t)
    ktot=kcloudy[0]+kcloudy[1] #Kabs + Kscat

    product= np.multiply(I_full,deltanu.as_matrix().ravel())
    norm=product.cumsum(axis=1).iloc[:,-1]
    
    kabs = np.multiply(product,kcloudy[0]).cumsum(axis=1).iloc[:,-1]
    kscat = np.multiply(product,kcloudy[1]).cumsum(axis=1).iloc[:,-1]
    prad= (pi/light_speed_micron)*np.multiply(product,ktot).cumsum(axis=1).iloc[:,-1]

    kabs=kabs/(norm) #erg cm-3
    kscat=kscat/(norm)
    
    prad=-prad
    kabs=(kabs.fillna(0))/(m_proton*10**n)
    kscat=(kscat.fillna(0))/(m_proton*10**n)
    
     
    ####
    
    file = open('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.kabs_I2'.format(directory,n,i,t),"w")
    file.write(kabs.to_string())
    file.close()
    
    file = open('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.kscat_I2'.format(directory,n,i,t),"w")
    file.write(kscat.to_string())
    file.close()
        
    file = open('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.prad'.format(directory,n,i,t),"w")
    file.write(prad.to_string())
    file.close()
     
    
    file = open('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.tau'.format(directory,n,i,t),"w")
    file.write(optical_depth.to_string())
    file.close()
    

    norm= (-norm)
    file = open('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.inc'.format(directory,n,i,t),"w")
    file.write(norm.to_string())
    file.close()


    return(prad,kabs,kscat,optical_depth,norm)
    
        
#%% 


for n in range(0,8):
    for i in logi:
        for t in logt:
                            print('doing n{0}_In{1:.1f}_Te{2:.2f}'.format(n,i,t))
                            kappa_I(n,i,t)

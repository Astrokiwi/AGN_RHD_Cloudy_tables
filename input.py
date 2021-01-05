#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:23:30 2017

@author: marta
"""

import numpy as np
import pandas as pd


directory='/Volumes/marta_filestore/cloudy_data/'
#directory='/home/marta/filestore/cloudy_data/'

#Intensity at the sublimation radius
I0=5.6*10**7

logi = []
rsub = []

#Set of intensities inside Rsub
logi.append(np.log10(I0*10**2.5))
logi.append(np.log10(I0*10**2))
logi.append(np.log10(I0*10**1.5))
logi.append(np.log10(I0*10))
logi.append(np.log10(I0*10**0.5))

#set of intensities further away Rsub
for i in range(14):
    a = I0 * 10**(-i*0.5)
    r = '{0:.2f}'.format(10**(i*0.5/2))
    logi.append(np.log10(a))
    rsub.append(r)

#Array of temperatures

logt=[]

T0=10
T_array=[]

for i in range(17):
    T=T0*10**(i*0.25)
    logt.append(np.log10(T))
    T_array.append(T)
    
    
for i in range(3):
    T = 10**(6+i)
    logt.append(np.log10(T))
    T_array.append(T)


#%%
#Taking a typical sample frequency set
n=0
i=logi[0]
t=8
frequencies=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.con'.format(directory,n,i,t),sep='\t',usecols=[0],names=['nu'],comment='#')

#%%
#Create a set of folders where the whole output from Cloudy will be stored
#==============================================================================
# import os
# for n in range(0,8):
#      for i in logi:
#        for t in logt:
#            #newpath='/Volumes/marta_filestore-1/cloudy_data/nodust/n{1}/In{2:.1f}/Te{3:.2f}/'.format(directory,n,i,t)
#            newpath='/Volumes/marta_filestore/cloudy_data/n{1}/In{2:.1f}/Te{3:.2f}/'.format(directory,n,i,t)

#            if not os.path.exists(newpath):
#                    os.makedirs(newpath)
#==============================================================================
#%%
def input_cloudy(I,n,temp):
    
    g = -2.8 + np.log10(1000)
    
    directory_file='{0}/n{1}/In{2:.1f}/Te{3:.2f}/'.format(directory,n,I,temp) 
   
    file = open("{0}n{2}_In{1:.1f}_Te{3:.2f}.in".format(directory_file,I,n,temp), "w") 
    
    string = '\
c === DEFINE CONTINUUM ===\n\
c\n\
AGN 5.5 -1.40 -0.50 -1.0\n\
Intensity {0:.2f} range total\n\
c\n\
c === fUV radiation field in 1000 G0 (Habing unit) === \n\
blackbody 30000 K\n\
intensity {1:.2f}, range 0.1 to 0.38 Ryd\n\
extinguish by 24, leakage = 0\n\
c\n\
c === (2) DEFINE PROPERTIES OF CLOUD ===\n\
c\n\
hden {2} \n\
abundances ism no grains\n\
grains ism function sublimation\n\
c\n\
CMB\n\
COSMIC RAYS BACKGROUND\n\
constant gas temperature {3} log \n\
stop column density 26\n\
iterate to convergence\n\
failures 2000 times\n\
stop temperature off\n\
c\n\
c === (3) CHANGE OUTPUT ===\n\
save cooling last ".cool"\n\
save dr last ".dr"\n\
save grain D/G ratio last ".ratio"\n\
save molecules last ".mol"\n\
save grain temperature last ".gtemp_full"\n\
save continuum last ".con"\n\
save pdr last ".pdr"\n\
save heating last ".het" \n\
'.format(I,g,n,temp)

    file.write(string)    

    for inu,nu in enumerate(frequencies.as_matrix()):
            file.write('save continuum emissivity {0} last ".znu{1}"\n\
'.format(nu[0],inu))
    
    file.close()
#%%

for n in range(0,8):
    for i in logi:
        for t in logt[-4:]:
            input_cloudy(i,n,t)



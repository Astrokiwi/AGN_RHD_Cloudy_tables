#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:28:01 2017

@author: marta
"""

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
import matplotlib as mpl
mpl.use('Agg')
import pylab as P



directory1="/Volumes/marta_filestore/cloudy_data/" #trillian
directory="/Volumes/marta_filestore-1/cloudy_data/" #somerville , I copied cool,pdr,prad2,tau_david,2,kabs2,kscat2,ratio,gtemp_full from srv1921
m_proton = 1.672621898*10**(-24)

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
#table = pd.read_csv("/Volumes/share/cooling_tables/tables_271117.txt",delimiter=r"\s+")
table = pd.read_csv(str(directory)+"070818.txt",delimiter=r"\s+")
#%%
table3 = pd.read_csv("/Volumes/share/cooling_tables/tables_271117.txt",delimiter=r"\s+")
#%%
logi= table.drop_duplicates('i0').i0.as_matrix()

logt=table.drop_duplicates('tgas').tgas.as_matrix()
colden=table.colden.iloc[:594]
#%%
table2=pd.DataFrame.copy(table)


#%%
def makedfnew(df,n,i,t):
    df2=df.query('tgas=={0} and i0=={1} and n0=={2}'.format(t,i,n))
    return(df2)
    #%%
def f(x,x_points,y_points):
    
     c=interpolate.interp1d(x_points, y_points, kind='nearest',fill_value='extrapolate')
     #c=interpolate.interp1d(x_points, y_points,fill_value='extrapolate')
     return c(x)
#%%
def f2(x,x_points,y_points,k):
    
     c=interpolate.interp1d(x_points, y_points, kind=k,fill_value='extrapolate')
     #c=interpolate.interp1d(x_points, y_points,fill_value='extrapolate')
     return c(x)
#%%
def queryneg(df):
    for col in df.columns[4:]:
         a=df[df[col]<0]
         if not a.empty:
             n1=a.drop_duplicates('n0').n0.as_matrix()
             i1= a.drop_duplicates('i0').i0.as_matrix()
             t1=a.drop_duplicates('tgas').tgas.as_matrix()
             print(col,n1,i1,t1)
             
#%%
for n in range(0,8):
    for i in logi: 
        for t in logt:
            df=makedfnew(table2,n,i,t)
            queryneg(df)
#%%
def querynan(df):
    for col in df.columns[4:]:
         a=df[col]
         if  a.isnull().values.any()==True:
            # n=a.drop_duplicates('n0').n0.as_matrix()
             #i= a.drop_duplicates('i0').i0.as_matrix()
            # t=a.drop_duplicates('tgas').tgas.as_matrix()
             print(col)

#table2.isin(['inf']).values.any()
#%%
#colden=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format("/Volumes/marta_filestore-1/cloudy_data/",2,logi[1],1),sep='\t',usecols=[1],names=['colden'],comment='#')  
#colden=colden.drop_duplicates('colden').reset_index(drop=True)
l=colden.shape[0]   
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
    #you need to normalize p_rad
   
   
    dgtot = dustToGasRatio.cumsum(axis=1).iloc[:,-1]
    dgtot=dgtot.as_matrix()
   
#==============================================================================
#     dg=df4.as_matrix()
#  
#     dgtot=[]
#      
#     for j in range(df4.shape[0]):
#          c=0       
#          for i in range(1,20):
#             c=c+dg[j,i]
#          dgtot.append(c)
#     dgtot=np.array(dgtot)
#==============================================================================
 
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
    df.loc[df['tau']<0,'tau']= f(df.query('tau<0').colden,df.query('tau>0').colden,df.query('tau>0').tau)# the first tau were negative
    return(df)

#%%
def lastcells(n,i,t,dfnew):
    
    old =cooling(n,i,t)
    old=pd.DataFrame(np.c_[np.full(old.shape[0],n),np.full(old.shape[0],i),np.full(old.shape[0],t),old.as_matrix()],columns=np.hstack((['n0','i0','tgas'],old.columns.values)))
    df=old.append(dfnew[dfnew.colden>old.colden.as_matrix()[-1]])  
    df=df.reset_index()

    tot=[]
    columns=np.array(['prad','kabs','kscat','tau'])

    for col in columns:
        c=interpolate.interp1d(df.colden,df[col],kind='linear')
        tot.append(c(colden))
    
    tot=np.c_[tot].T[0]


    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
  
    return(d)
#%%   
def df_interpolated_smoothed(df,n,i,t):
    
    columns = df.columns[4:]
    
    h2=f(colden,df.colden,df.heat)
    c2=f(colden,df.colden,df.cool)
    tg2=f(colden,df.colden,df.tgrain)
    p=f(colden,df.colden,df.arad)
    dg=f(colden,df.colden,df.dg)
    tau=f(colden,df.colden,df.tau)

    
    h2=signal.savgol_filter(h2,11,1,mode='nearest')
    c2=signal.savgol_filter(c2,11,1,mode='nearest')
    tg2=signal.savgol_filter(tg2,11,1,mode='nearest')
    p=signal.savgol_filter(p,11,1,mode='nearest')
    dg=signal.savgol_filter(dg,9,1,mode='nearest')
    tau=signal.savgol_filter(tau,11,1,mode='nearest')

    
    
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,h2,c2,tg2,p,dg,tau],columns=['n0','i0','tgas','colden','heat','cool','tgrain','arad','dg','tau'])
    
    return(d)
    #%%

def InterpForT(n,i,t): # t_wanted is the index 
    
    t_wanted=np.where(logt==t)[0][0]
    a=[]
    for ti in logt[0:t_wanted]:
            df=makedfnew(table2,n,i,ti)
            a.append(df)
     
    d=pd.DataFrame(np.vstack(a),columns=df.columns.values)

    tot2=[]

    columns=np.array(['prad','kabs','kscat','tau'])
    for col in columns:
        tot=[]
        for N in colden.as_matrix():
            df=d.query('colden=={0}'.format(N))
            tot.append(f(t,df.tgas,df[col])) 

        tot2.append(tot)
        
    tot2=np.c_[tot2].T    
    
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
    
    return(d)
#%%
def InterpForT_reversed(n,i,t):
    
    logt_rev=logt[::-1]
    t_wanted=np.where(logt_rev==t)[0][0]
    
    a=[]
    
    for ti in logt_rev[3:t_wanted]: #array([ 5.  ,  4.75,  4.5 ,  4.25,  4.  ,  3.75,  3.5 ...])
            df=makedfnew(table2,n,i,ti)
            a.append(df)
       
    d=pd.DataFrame(np.vstack(a),columns=df.columns.values)

    tot2=[]
    columns=np.array(['prad','kabs','kscat','tau'])
    for col in columns:
        tot=[]
        for N in colden.as_matrix():
            df=d.query('colden=={0}'.format(N))
            tot.append(f(t,df.tgas,df[col])) 

        tot2.append(tot)
        
    tot2=np.c_[tot2].T    
    
   # d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=df.columns.values)
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
    
#==============================================================================
#     -
#==============================================================================
    
    return(d)
#%%
def ExtrapByn(n,i,t):
    
    a=[]
    for ni in range(0,n):
            df= makedfnew(table2,ni,i,t)
            a.append(df)
       
    d=pd.DataFrame(np.vstack(a),columns=df.columns.values)

    tot2=[]
    columns=['kabs','kscat','prad','tau']
    #columns=['tau']
    for col in columns:
        tot=[]
        for N in colden.as_matrix():
            df=d.query('colden=={0}'.format(N))
            tot.append(f(10**n,10**df.n0,df[col])) 

        tot2.append(tot)
        
    tot2=np.c_[tot2].T    

    
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
    
    return(d) 
#%%
def ExtrapByI(n,i,t):
    l=len(colden)
    i_wanted=np.where(logi==i)[0][0]
    a=[]
    for ii in logi[:i_wanted]:
            df=makedfnew(table2,n,ii,t)
            a.append(df)
            
    if (n==7)and(i==logi[-2]):
        a.append(makedfnew(table2,n,logi[-1],t))
    #if (n==7)and(i==logi[-3]):
    #    a.append(makedfnew(table2,n,logi[-2],t))
     #   a.append(makedfnew(table2,n,logi[-1],t)) 
    
    d=pd.DataFrame(np.vstack(a),columns=['n0','i0','tgas','colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
   
    tot2=[]
    columns=['kabs','kscat','prad','tau']

    for col in columns:
    #for col in np.array(['tau']):
        tot=[]
        for N in colden.as_matrix():
            df=d.query('colden=={0}'.format(N))
            df=df.sort_values('i0', axis=0, ascending=True)
            tot.append(f(i,df.i0,df[col])) 

        tot2.append(tot)
        
    tot2=np.c_[tot2].T    
  
    #d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=['n0','i0','tgas','colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot2],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
    return(d)
#%%
n=4
#==============================================================================
# 4 2.248188 2.0 6.29e+25   interp for T
# 4 2.248188 2.25 3.41e+25 interp for T
# 4 1.748188 2.0 8.02e+25 interp for T
# 4 1.748188 2.25 5.67e+25 interp for T
# 4 1.248188 2.0 7.37e+25
# 4 1.248188 2.25 7.27e+25
#==============================================================================
for i in logi[-3:]:
    for t in logt[::-1]:
        a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
        if a[0]<10**26:
               # print(n,i,t,a[0])
                df=InterpForT_reversed(n,i,t)
                #df=makedfnew(table,n,i,t)
                mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
                table2.loc[mask,'tau']=df.tau.as_matrix().ravel()
                
                df2=cooling(n,i,t)

                oldvsnew(df,df2)
#%%
t=2.25
i=logi[-3]
df=ExtrapByn(n,i,t)

mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
table2.loc[mask,'tau']=df.tau.as_matrix().ravel()
#%%
def oldvsnew(df,df2):
                i=df.drop_duplicates('i0').i0.as_matrix()[0]
                t=df.drop_duplicates('tgas').tgas.as_matrix()[0]
                n=df.drop_duplicates('n0').n0.as_matrix()[0]
    
                fig = plt.figure()
                fig.suptitle("n{0}_I{1:.1f}_Te{2}".format(n,i,t), color='green',fontsize=16,fontweight='bold')
                   
                ax = plt.subplot("221")
                ax.set_title('non interp tau')
                ax.semilogx(df2.colden,df2.tau,'blue')
                
                ax = plt.subplot("222")
                ax.set_title('tau')
                ax.semilogx(df.colden,df.tau,'red')
                 
                ax = plt.subplot("223")
                ax.set_title(' non interp prad')
                ax.semilogx(df2.colden,df2.prad,'blue')
                
                ax = plt.subplot("224")
                ax.set_title('prad')
                ax.semilogx(df.colden,df.prad,'red')
                plt.show()
                                              
#%%
#==============================================================================
# 5 4.748188 2.75 6.62e+25 #t
# 5 4.248188 2.75 4.56e+25 #t
# 5 3.748188 2.75 2.01e+25 #t
# 5 3.248188 2.0 7.71e+25 #n
# 5 3.248188 2.25 4.81e+25 #
# 5 3.248188 2.75 8.19e+24 #t
# 5 2.748188 1.0 8.34e+25
# 5 2.748188 1.25 7.24e+25
# 5 2.748188 1.5 8.65e+25
# 5 2.748188 1.75 3.27e+25
# 5 2.748188 2.0 3.13e+25 #n
# 5 2.748188 2.25 1.63e+25 #t
# 5 2.748188 2.75 2.79e+24 #t
# 5 2.248188 1.0 5.42e+25
# 5 2.248188 1.25 4.62e+25
# 5 2.248188 1.5 3.25e+25
# 5 2.248188 1.75 2.28e+25   INTERPOLA PER T AL CONTRARIO (in totale sono 36)
# 5 2.248188 2.0 1.43e+25 #n
# 5 2.248188 2.25 9.42e+24
# 5 2.248188 2.5 3.56e+25  #i
# 5 2.248188 2.75 3.91e+23 #t
# 5 1.748188 1.0 5.04e+25
# 5 1.748188 1.25 3.9e+25
# 5 1.748188 1.5 3.03e+25
# 5 1.748188 1.75 2.2e+25
# 5 1.748188 2.0 1.84e+25 #n
# 5 1.748188 2.25 1.26e+25 #t
# 5 1.748188 2.75 9.11e+23 #t
# 5 1.248188 1.0 4.67e+25
# 5 1.248188 1.25 3.89e+25
# 5 1.248188 1.5 3.79e+25
# 5 1.248188 1.75 3.6e+25
# 5 1.248188 2.0 3.9e+25 #n
# 5 1.248188 2.25 2.96e+25 #t
# 5 1.248188 2.75 2.38e+24 #t
# 5 1.248188 3.0 9.84e+25 #no doing it
#==============================================================================

n=5
t=2.75

#mask=df.colden>10**22
for i in logi:
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                print(n,i,t,a[0])
                df=InterpForT_reversed(n,i,t)
                df2=cooling(n,i,t)
                
                oldvsnew(df,df2)
                mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
                
                
                for col in np.array(['prad','kabs','kscat','tau']):
                    table2.loc[mask,col]=df[col].as_matrix().ravel()
                
#%%
n=5
t=2.25
for i in logi:
    if i!=logi[-3]:
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
 
                    print(n,i,t,a[0])
                    df=InterpForT_reversed(n,i,t) #5.  ,  4.75,  4.5 ,  4.25,  4.  ,  3.75,  3.5 ,  3.25,  3.  , 2.75,  2.5 ]
                    df2=cooling(n,i,t)
                    oldvsnew(df,df2)
                    mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
                    for col in np.array(['prad','kabs','kscat','tau']):
                        table2.loc[mask,col]=df[col].as_matrix().ravel()
                   

#%%
n=5
i=logi[-3]
t=2.25
df=ExtrapByI(n,i,2.25)

oldvsnew(df,df2)
mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
for col in np.array(['prad','kabs','kscat','tau']):
                        table2.loc[mask,col]=df[col].as_matrix().ravel()
                        
#%%
n=5

for i in logi: 
    for t in logt[0:5]: #logt=1,1.25,1.75
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:    
                    print(n,i,t,a[0])
                    df2=cooling(n,i,t)
                    df=ExtrapByn(n,i,t)
                    mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n)
                    table2.loc[mask,'tau']=df.tau.as_matrix().ravel()
                    
                    fig = plt.figure()
                    fig.suptitle("n{0}_I{1:.1f}_Te{2}".format(n,i,t), color='green',fontsize=16,fontweight='bold')
                       
                    ax = plt.subplot("121")
                    ax.set_title('non interp tau')
                    ax.semilogx(df2.colden,df2.tau,'blue')
                    
                    ax = plt.subplot("122")
                    ax.set_title('tau')
                    ax.semilogx(df.colden,df.tau,'red')
                    plt.show()
#%%

c=0
n=6
#mask=df.colden>10**22
for i in logi:
        for t in logt:
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 c=c+1
                 print(n,i,t,a[0])
#==============================================================================
# 6 4.248188 2.25 9.35e+25
# 6 3.748188 1.25 9.12e+25
# 6 3.748188 1.5 6.15e+25
# 6 3.748188 1.75 7.69e+25
# 6 3.748188 2.0 3.89e+25
# 6 3.748188 2.25 2.69e+25
# 6 3.248188 1.0 6.22e+25
# 6 3.248188 1.25 4.49e+25
# 6 3.248188 1.5 3.57e+25
# 6 3.248188 1.75 3.96e+25
# 6 3.248188 2.0 2.77e+25
# 6 3.248188 2.25 1.63e+25
# 6 2.748188 1.0 1.42e+25
# 6 2.748188 1.25 2.86e+24
# 6 2.748188 1.5 2.64e+24
# 6 2.748188 1.75 1.92e+24
# 6 2.748188 2.0 1.65e+24
# 6 2.748188 2.25 2.21e+24
# 6 2.748188 2.5 7.28e+25
# 6 2.248188 1.0 1.36e+24
# 6 2.248188 1.25 8.81e+23
# 6 2.248188 1.5 7.77e+23
# 6 2.248188 1.75 7.29e+23
# 6 2.248188 2.0 5.72e+23
# 6 2.248188 2.25 4.88e+23
# 6 2.248188 2.5 1.51e+25
# 6 1.748188 1.0 8.12e+24
# 6 1.748188 1.25 7.11e+24
# 6 1.748188 1.5 8.26e+24
# 6 1.748188 1.75 4.63e+24
# 6 1.748188 2.0 1.9e+25
# 6 1.748188 2.25 8.75e+24
# 6 1.748188 3.0 7.58e+25
# 6 1.248188 1.0 5.74e+25
# 6 1.248188 1.25 5.78e+25
# 6 1.248188 1.5 4.64e+25
# 6 1.248188 1.75 3.94e+25
# 6 1.248188 2.0 3.34e+25
# 6 1.248188 2.25 1.68e+25
# 6 1.248188 3.0 6.99e+25
#==============================================================================
#%%
n=6

#mask=df.colden>10**22
for t in np.array([3,2.5,2.25]):
    for i in logi:
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 print(n,i,t,a[0])
                 df=InterpForT_reversed(n,i,t)
                 df2=cooling(n,i,t)
                 oldvsnew(df,df2)
                 mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                 for col in np.array(['prad','kabs','kscat','tau']):
                        table2.loc[mask,col]=df[col].as_matrix().ravel()
#%%6 2.248188 2.25 4.88e+23
n=6
i=logi[-3]
t=2.25
df=ExtrapByI(n,i,t)

mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
for col in np.array(['prad','kabs','kscat','tau']):
                        table2.loc[mask,col]=df[col].as_matrix().ravel()
#%%
c=0
n=6
#mask=df.colden>10**22
for i in logi:
    for t in np.array([2,1.75,1.5,1.25,1]):
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 c=c+1
                 print(n,i,t,a[0])
                 df=ExtrapByn(n,i,t)
                 df2=cooling(n,i,t)
                 oldvsnew(df,df2)
                 mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                 for col in np.array(['prad','kabs','kscat','tau']):
                       table2.loc[mask,col]=df[col].as_matrix().ravel()
                     

#%%   
'''
n=6        
i=logi[-3]
t=1.5

df=InterpForT(n,i,t)
mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
for col in np.array(['prad','kabs','kscat','tau']):
                       table2.loc[mask,col]=df[col].as_matrix().ravel()

#%%            
n=6
#mask=df.colden>10**22
for i in logi[-4:-1]:
    for t in np.array([1.75,2]):
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 print(n,i,t,a[0])
                 df=InterpForT(n,i,t)
                 df2=cooling(n,i,t)
                 oldvsnew(df,df2)
#%%
c=0
n=6
#mask=df.colden>10**22
i=logi[-5]
for t in np.array([2,1.75,1.5,1.25,1]):
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 c=c+1
                 print(n,i,t,a[0])
                 df=ExtrapByI(n,i,t)
                 df2=cooling(n,i,t)
                 oldvsnew(df,df2)
                 mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                 for col in np.array(['prad','kabs','kscat','tau']):
                      table2.loc[mask,col]=df[col].as_matrix().ravel()
'''
#%%
def fixlastcells(n,i,t):
    df=cooling(n,i,t)
    #dfnew=makedfnew(table,n,i,t)
    df_last=df.iloc[-30:,:]
    colden2=colden[colden.colden>df_last.colden.as_matrix()[-1]]
    coldentot=np.concatenate((df.colden.as_matrix().ravel(),colden2.as_matrix().ravel()),axis=0)
    
    tot=[]

    for col in df.columns[1:]:
         a=df[col].as_matrix()
         a=np.concatenate((a,f2(colden2,df_last.colden,df_last[col],'nearest').ravel()))
         
        # cs=interpolate.InterpolatedUnivariateSpline(df_last.colden,df_last[col],k='1')
         #a=np.concatenate((a,cs(colden2).ravel()),axis=0)
        # cs=interpolate.InterpolatedUnivariateSpline(coldentot,a,k='1')
         
         #tot.append(cs(colden))
         tot.append(f(colden,coldentot,a))
    
    tot=np.c_[tot].T[0]

    
    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot],columns=['n0','i0','tgas','colden','tgrain','heat','cool','prad','dg','kabs','kscat','tau'])
    #d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,dfnew[df.columns[1:-4]],tot],columns=np.hstack((['n0','i0','tgas'],df.columns.values)))

    return(d)
#%%
n=7
c=0

#mask=df.colden>10**22
for i in logi:
    for t in logt:
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                 c=c+1
                 print(n,i,t,a[0])
                 
#==============================================================================
# 7 6.248188 2.75 9.56e+25
# 7 4.248188 1.0 5.49e+25
# 7 4.248188 1.25 4.59e+25
# 7 4.248188 1.5 4.26e+25
# 7 4.248188 1.75 5.76e+25
# 7 4.248188 2.0 7.27e+25
# 7 3.748188 1.0 4.19e+25
# 7 3.748188 1.25 2.54e+25
# 7 3.748188 1.5 2.35e+25
# 7 3.748188 1.75 2.44e+25
# 7 3.748188 2.0 4.16e+25
# 7 3.748188 2.25 8.05e+25
# 7 3.248188 1.0 3.44e+25
# 7 3.248188 1.25 2.71e+25
# 7 3.248188 1.5 1.95e+25
# 7 3.248188 1.75 1.42e+25
# 7 3.248188 2.0 3.34e+25
# 7 3.248188 2.25 5.83e+25
# 7 2.748188 1.0 5.51e+24
# 7 2.748188 1.25 2.76e+24
# 7 2.748188 1.5 2.32e+24
# 7 2.748188 1.75 3.17e+24
# 7 2.748188 2.0 5.42e+24
# 7 2.748188 2.25 5.02e+25
# 7 2.748188 3.0 7.88e+25
# 7 2.248188 1.0 4.55e+24
# 7 2.248188 1.25 3.23e+24
# 7 2.248188 1.5 3.88e+24
# 7 2.248188 1.75 4.58e+24
# 7 2.248188 2.0 1.13e+25
# 7 2.248188 2.25 5.47e+25
# 7 1.748188 1.0 1e+25
# 7 1.748188 1.25 1.23e+25
# 7 1.748188 1.5 3.96e+24
# 7 1.748188 1.75 8.98e+24
# 7 1.748188 2.0 8.71e+25
# 7 1.248188 1.0 3.83e+25
# 7 1.248188 1.25 3.16e+25
# 7 1.248188 1.5 1.5e+24
# 7 1.248188 1.75 4.56e+25
#==============================================================================
#%%
n=7
for i in logi:
    for t in np.array([3,2.75,2.25,2]):
            a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
            if a<10**26:
                print(n,i,t,a[0])

                df2=cooling(n,i,t)
                df=InterpForT_reversed(n,i,t)
                oldvsnew(df,df2)
                mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                for col in np.array(['prad','kabs','kscat','tau']):
                      table2.loc[mask,col]=df[col].as_matrix().ravel()

#%%
n=7
t=1.75
for i in logi:
   # for t in np.array([1,1.25,1.5,1.75]):
       # if i=!logi[-1]
        a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
        if a<10**26:
                print(n,i,t,a[0])
                df2=cooling(n,i,t)
                df=InterpForT_reversed(n,i,t)
                #df=lastcells(n,i,t,df)
                oldvsnew(df,df2)
                plt.loglog(colden,df.prad)
                mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                for col in np.array(['prad','kabs','kscat','tau']):
                   table2.loc[mask,col]=df[col].as_matrix().ravel()
      #  else:
       #        df=ExtrapByI(n,i,t)
               
      #  mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
       # for col in np.array(['prad','kabs','kscat','tau']):
        #           table2.loc[mask,col]=df[col].as_matrix().ravel()
#%%
n=7
t=1.75
for i in logi:
    df=makedfnew(table2,n,i,t)
    plt.loglog(df.colden,df.prad,label='{0:.1f}'.format(i))
plt.show()
#%%
n=7
#1.5 fatto
for t in logt[0:2]:
    for i in logi[1:]:
     a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
     if a<10**26:
        print(n,i,t,a[0])
        df=ExtrapByI(n,i,t)
        mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
        for col in np.array(['prad','kabs','kscat','tau']):
            cs=interpolate.pchip(colden.as_matrix().ravel(),df[col])
            table2.loc[mask,col]=cs(colden).ravel()
        #Edf=lastcells(n,i,t,df)
        df=makedfnew(table2,n,i,t)
        df2=cooling(n,i,t)
        oldvsnew(df,df2)
        plt.loglog(df.colden,df.prad,label='{0:.1f}'.format(i))
        plt.show()
#%%
n=7
t=1.25
for i in logi[1:]:
    #if (i!=logi[-1]):
        a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
        if a<10**26:
                print(n,i,t,a[0])
                df=InterpForT_reversed(n,i,t)
                #df=lastcells(n,i,t,df)
                df2=cooling(n,i,t)
                oldvsnew(df,df2)
                plt.loglog(df.colden,df.prad,label='{0:.1f}'.format(i))
                plt.show()
#%%
n=7
for t in np.array([1.5,1.25]):
    for i in logi[1:]: 
     a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
     if a<10**26:
        df=InterpForT_reversed(n,i,t)
       # df=lastcells(n,i,t,df)
        mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
        for col in np.array(['prad','kabs','kscat','tau']):
                   table2.loc[mask,col]=df[col].as_matrix().ravel()


#%%
n=7
t=1
for i in np.array([logi[-5],logi[-3]]):    
     a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
     if a<10**26:
          print(n,i,t,a[0])
          df=ExtrapByI(n,i,t)
                #df=lastcells(n,i,t,df)
          df2=cooling(n,i,t)
          oldvsnew(df,df2)
          mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
          for col in np.array(['prad','kabs','kscat','tau']):
                   table2.loc[mask,col]=df[col].as_matrix().ravel()

#%%
t=1
i=logi[-3]
df=ExtrapByn(n,i,t)
df2=cooling(n,i,t)
oldvsnew(df,df2)

#%%
t=1
i=logi[-1]
df=ExtrapByn(n,i,t)
#df2=cooling(n,i,t)
#oldvsnew(df,df2)
#%%

   
mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
for col in np.array(['prad','kabs','kscat','tau']):
                   table2.loc[mask,col]=df[col].as_matrix().ravel()
#plt.legend()

#%%
n=7
#columns=['heat','cool','prad','kabs','kscat','tau']
columns=['tgrain','dg']

for n in range(0,8):
    for i in logi:
        for t in logt:
            df=makedfnew(table,n,i,t).drop_duplicates('colden')
            mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
            for col in columns:
                       table2.loc[mask,col]=df[col].as_matrix().ravel()

#%%
#TO FIX THE NEGATIVE VALUES FOUND with the function queryneg

for n in np.array([0,1,7,]):
    for i in logi[(logi>7.7)]:
        for t in logt:
            df=makedfnew(table2,n,i,t)
            mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
            
            for col in df.columns[4:]:
                a=df[df[col]<0]
                if not a.empty:
                    df.loc[df[col]<0,col]= f(df.query(col+'<0').colden,df.query(col+'>0').colden,df.query(col+'>0')[col])
                    table2.loc[mask,col]=df[col].as_matrix().ravel()
                    
                    
''''
i=logi[-2]
for t in np.array([1,1.25,1.5,1.75]):
        a= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[1],names=['colden'],comment='#').as_matrix()[-1]
        if a<10**26:
                print(n,i,t,a[0])
                df2=cooling(n,i,t)
                df=ExtrapByI(n,i,t) #EXTRAPOLATEBYI E' STATA MODIFICATA AGGIUNGENDO LOGI[-1]
                df3=ExtrapByn(n,i,t)
                print('OLD_VS_NEW_BY_I')
                oldvsnew(df,df2)
                print('OLD_VS_NEW_BY_N')
                oldvsnew(df3,df2)
                
                mask=(table2.i0==i)&(table2.tgas==t)&(table2.n0==n) 
                table2.loc[mask,'tau']=df['tau'].as_matrix().ravel()
                for col in np.array(['prad','kabs','kscat']):
                      table2.loc[mask,col]=df3[col].as_matrix().ravel()
    
'''
#%%left with
#==============================================================================
# 7 4.248188 1.0 5.49e+25
# 7 4.248188 1.25 4.59e+25
# 7 4.248188 1.5 4.26e+25
# 7 4.248188 1.75 5.76e+25
# 7 3.748188 1.0 4.19e+25
# 7 3.748188 1.25 2.54e+25
# 7 3.748188 1.5 2.35e+25
# 7 3.748188 1.75 2.44e+25
# 7 3.248188 1.0 3.44e+25
# 7 3.248188 1.25 2.71e+25
# 7 3.248188 1.5 1.95e+25
# 7 3.248188 1.75 1.42e+25
# 7 2.748188 1.0 5.51e+24
# 7 2.748188 1.25 2.76e+24
# 7 2.748188 1.5 2.32e+24
# 7 2.748188 1.75 3.17e+24
# 7 2.248188 1.0 4.55e+24
# 7 2.248188 1.25 3.23e+24
# 7 2.248188 1.5 3.88e+24
# 7 2.248188 1.75 4.58e+24
#==============================================================================
    
file = open('{0}tables_070818.txt'.format('/Volumes/share/cooling_tables/'),"w")
#file= open('{0}131117.txt'.format(directory),"w")
file.write(table2.to_string())
file.close()#%

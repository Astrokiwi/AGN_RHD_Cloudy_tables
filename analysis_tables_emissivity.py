#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 16:41:43 2020

@author: marta
"""

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
import matplotlib as mpl
import pylab as P



directory="/Volumes/marta_filestore-1/cloudy_data/" #somerville , I copied cool,pdr,prad2,tau_david,2,kabs2,kscat2,ratio,gtemp_full from srv1921


#%%
#table = pd.read_csv("/Volumes/share/cooling_tables/tables_271117.txt",delimiter=r"\s+")
table = pd.read_csv('/Users/marta/Desktop/'+"emissivity.txt",delimiter=r"\s+")

#%%
logi= table.drop_duplicates('i0').i0.as_matrix()

logt=table.drop_duplicates('tgas').tgas.as_matrix()
#%%
colden=table.colden.iloc[:594]
#%%
table3=pd.DataFrame.copy(table)


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
            df=makedfnew(table3,n,i,t)
            queryneg(df)

#table2.isin(['inf']).values.any()
#%%
#colden=pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format("/Volumes/marta_filestore-1/cloudy_data/",2,logi[1],1),sep='\t',usecols=[1],names=['colden'],comment='#')  
#colden=colden.drop_duplicates('colden').reset_index(drop=True)
l=colden.shape[0]   
#%%
def emissivity(n,i,t):

    colden= pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.pdr'.format(directory,n,i,t),sep='\t',usecols=[0,1],names=['depth','colden'],comment='#')
    
    nujnu12 = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,3616),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values
    nujnu18 =  pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,3697),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values
    nujnu850 = pd.read_csv('{0}n{1}/In{2:.1f}/Te{3:.2f}/n{1}_In{2:.1f}_Te{3:.2f}.znu{4}'.format(directory,n,i,t,4468),sep='\t',usecols=[3],names=['nujnu'],skiprows=[0],comment='#').values

    table = np.c_[colden.colden,nujnu12,nujnu18,nujnu850]
        #table=np.c_[colden.colden,cool.Htot,cool.Cool,prad,kabs_I,kscat_I,tau_david]

    df=pd.DataFrame(table,columns=['colden','nujnu12','nujnu18','nujnu850'])
    #df=pd.DataFrame(table,columns=['colden','heat','cool','prad','kabs','kscat','tau'])
    return(df)
#%%
def lastcells(n,i,t,dfnew):
    
    old =emissivity(n,i,t)
    old=pd.DataFrame(np.c_[np.full(old.shape[0],n),np.full(old.shape[0],i),np.full(old.shape[0],t),old.as_matrix()],columns=np.hstack((['n0','i0','tgas'],old.columns.values)))
    df=old.append(dfnew[dfnew.colden>old.colden.as_matrix()[-1]])  
    df=df.reset_index()

    tot=[]
    columns=np.array(['nujnu12','nujnu18','nujnu850'])

    for col in columns:
        c=interpolate.interp1d(df.colden,df[col],kind='linear')
        tot.append(c(colden))
    
    tot=np.c_[tot].T[0]


    d=pd.DataFrame(np.c_[np.full(l,n),np.full(l,i),np.full(l,t),colden,tot],columns=np.hstack((['n0','i0','tgas','colden'],columns)))
        
    return(d)

    #%%

def InterpForT(n,i,t): # t_wanted is the index 
    
    t_wanted=np.where(logt==t)[0][0]
    a=[]
    for ti in logt[0:t_wanted]:
            df=makedfnew(table3,n,i,ti)
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
            df=makedfnew(table3,n,i,ti)
            a.append(df)
       
    d=pd.DataFrame(np.vstack(a),columns=df.columns.values)

    tot2=[]
    columns=np.array(['nujnu12','nujnu18','nujnu850'])
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
            df= makedfnew(table3,ni,i,t)
            a.append(df)
       
    d=pd.DataFrame(np.vstack(a),columns=df.columns.values)

    tot2=[]
    columns=['nujnu12','nujnu18','nujnu850']
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
            df=makedfnew(table3,n,ii,t)
            a.append(df)
            
    if (n==7)and(i==logi[-2]):
        a.append(makedfnew(table3,n,logi[-1],t))
    #if (n==7)and(i==logi[-3]):
    #    a.append(makedfnew(table2,n,logi[-2],t))
     #   a.append(makedfnew(table2,n,logi[-1],t)) 
    
    d=pd.DataFrame(np.vstack(a),columns=['n0','i0','tgas','colden','nujnu12','nujnu18','nujnu850'])
   
    tot2=[]
    columns=['nujnu12','nujnu18','nujnu850']

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
def oldvsnew(df,df2):
                i=df.drop_duplicates('i0').i0.as_matrix()[0]
                t=df.drop_duplicates('tgas').tgas.as_matrix()[0]
                n=df.drop_duplicates('n0').n0.as_matrix()[0]
    
                fig = plt.figure()
                fig.suptitle("n{0}_I{1:.1f}_Te{2}".format(n,i,t), color='green',fontsize=16,fontweight='bold')
                   
                ax = plt.subplot("221")
                ax.set_title('non interp 12')
                ax.semilogx(df2.colden,df2.nujnu12,'blue')
                
                ax = plt.subplot("222")
                ax.set_title('12')
                ax.semilogx(df.colden,df.nujnu12,'red')
                 
                ax = plt.subplot("223")
                ax.set_title(' non interp 850')
                ax.semilogx(df2.colden,df2.nujnu850,'blue')
                
                ax = plt.subplot("224")
                ax.set_title('nujnu850')
                ax.semilogx(df.colden,df.nujnu850,'red')
                plt.show()
                                              


#
#%%

n=7
t=1.25
for i in logi:
    df= makedfnew(table3,n,i,t)
    plt.loglog(df.colden,df.nujnu18,label=i)
plt.legend()
plt.close('All')
#plt.legend()
#%%
i=logi[-7]
df2 = InterpForT_reversed(n,i,t)
#%%

mask=(table3.i0==i)&(table3.tgas==t)&(table3.n0==n) 
for col in np.array(['nujnu18']):
                   table3.loc[mask,col]=df2[col].as_matrix().ravel()

  #%%  
file = open('{0}emissivity_140620.txt'.format('/Users/marta/Desktop/'),"w")
#file= open('{0}131117.txt'.format(directory),"w")
file.write(table3.to_string())
file.close()#%
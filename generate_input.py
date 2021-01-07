#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:23:30 2017

@author: marta
"""

import numpy as np
import os
import parameter_ranges as pr

def input_cloudy_simple(I,n,temp,frequencies,template):

    g = -2.8 + np.log10(1000)

    file_dir=f'{pr.directory}/n{n:.0f}/In{I:.1f}/Te{temp:.2f}/'
    file_name=f"{file_dir}n{n:.0f}_In{I:.1f}_Te{temp:.2f}.in"

    with open(file_name, "w") as file:
        file.write(
            template.format(I=I, n=n, temp=temp, g=g)
        )

        for inu,nu in enumerate(frequencies):
            file.write(f'save continuum emissivity {nu} last units microns ".znu{inu}"\n')

#%%

if __name__=='__main__':

    # %%
    # Create a set of folders where the whole output from Cloudy will be stored
    # ==============================================================================
    for n in pr.logn :
        for i in pr.logi :
            for t in pr.logt :
                newpath = f"{pr.directory}/n{n:.0f}/In{i:.1f}/Te{t:.2f}"

                if not os.path.exists(newpath) :
                    os.makedirs(newpath)
    # ==============================================================================
    # %%

    #%%
    #Taking a typical sample frequency set
    frequencies=np.loadtxt('continuum_freq.txt')

    with open("input_template.txt") as f:
        template = f.read()

    for n in pr.logn[0:2]:
        for i in pr.logi[0:2]:
            for t in pr.logt[0:2]:
                print(n,i,t)
                input_cloudy_simple(i,n,t,frequencies,template)



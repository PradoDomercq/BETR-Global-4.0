################################################################################
'''Module "mkDmixQ", JUNE 2017, is part of
BETR-Research by Fangyuan Zhao

This module calculates Dmix values with Q for mixing processes of a
particular model parametrization.'''
################################################################################
import numpy
from numpy import *
import inspect
from globalz import *
from helpers import *
import sys
import copy
import scipy.sparse as sp
import re
import pdb
import write_ncfile
import os

def mkDmixQ(m):      
    comp=m.nocomp;     cells=m.nocells;  seasons=m.nots
    vz=zeros((comp,cells,seasons));  T=zeros((2,cells,seasons))
    dk=zeros((comp,cells,seasons));  Kloss=zeros((2,cells,seasons))
    Qmix=zeros((2,cells,seasons))

    ''' uses zdict and vdict to construct vz'''
    for c in m.compdict.keys():
        vz[c-1]=m.zdict[c]['bulk']*m.vdict[c]['bulk']
        
    '''T in mixing process'''        
    for [k,v] in m.Dproc.items():
        key=k[0]-1
        mstring = k[2]
        if re.search("mix", mstring): 
           if key==0: T[0]=m.par['h1']/m.par['mixing12']
           if key==1: T[1]=m.par['h2']/m.par['mixing21']

    '''Kloss'''  
    for [k,v] in m.Dproc.items():
        key=k[0]-1
        mstring = k[2]
        if re.search("deg|wet|disso", mstring): 
           if key==0: dk[0]+=v
           if key==1: dk[1]+=v
           
    for i in range(2):       
        for ii in range(cells):
            for iii in range(seasons):
                Kloss[i][ii][iii]=dk[i][ii][iii]/vz[i][ii][iii]
                
    Qmix=(1+Kloss*T)/exp(Kloss*T)


    return(Qmix)
 
#################################################################################
'''  
                
    Qmix=Kloss*T    
    timesteps=12
    if m.nocells==288:
       nlat=12;   nlon=24
    if m.nocells==1152:
       nlat=24;   nlon=48
    if m.nocells==4608:
       nlat=48;   nlon=96
    outKTeast=zeros((timesteps,2,nlat,nlon))
    fn1=os.path.join('Output', m.run, 'hr_KTmix_24000.nc')   
    for c in range(2):
        for t in range(timesteps):
            outKTeast[t,c,:,:]=Qmix[c][:,t].reshape(nlat,nlon)            
    write_ncfile.writenc(outKTeast, fn1, True, True, 1, varname='V',unit='mol')
    
    pdb.set_trace()
'''

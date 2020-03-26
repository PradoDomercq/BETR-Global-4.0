# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 17:42:28 2013
LinFluxApprox provides helpers needed for the Linear Approximation of 
calculated inter and intramedia Fluxes implemented in BETRS 3.1

@author: SSchenker
"""
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
from mkbigD import *      #FY

def mkFluxkey(m, flowmat, procmat, lat=12, lon=24, use_correction=False):
    keylist=[]
    for season in range(12):
        keylist.append([])
        for comp in range(lat * lon):
            keylist[season].append([])
            for medium in range(10):
                keylist[season][comp].append({})

    for pkey in procmat.keys():
        key, medium = match_keys(pkey)
        for comp in range (len(procmat[pkey])):
            for season in range(len(procmat[pkey][comp])):                
                if key in keylist[season][comp][medium].keys():
                    keylist[season][comp][medium][key] += procmat[pkey][comp][season]
                else:
                    keylist[season][comp][medium][key] = procmat[pkey][comp][season]

    if use_correction==True:
        '''The following calculate 'pseudo dilution'=D-DQn+1 and flow_Q transport which instead the flow transport. added by Fangyuan, 2rd Dec., 2016'''                    
        flowmatQ=mkflowQ(m)
        #calculate flow_Q transport 
        for fkeyQ in  flowmatQ.keys():
            key="flux"
            medium = fkeyQ[0]-1
            for entry in flowmatQ[fkeyQ]:
                comp = int(entry[0])-1            
                for i in range (2,14):
                    season=i-2
                    if key in keylist[season][comp][medium].keys():
                        keylist[season][comp][medium][(key,entry[1])] += entry[i]
                    else:
                        keylist[season][comp][medium][(key,entry[1])] = entry[i]                    
                    
        # remove intra-cell flows
        for [k,v] in flowmat.items():
            intercell=list(where(v[:,0]!=v[:,1])[0])
            intracell=list(where(v[:,0]==v[:,1])[0])
            flowmat[k]=v[intercell,:]    
        #calculate 'pseudo dilution'=D-DQn+1
        for [k,v] in flowmatQ.items():
            for [kk,vv] in flowmat.items():
                if k==kk:
                    v[:,2:14]=vv[:,2:14]-v[:,2:14]   
        #change pseudo dilution in keylist form
        for fkeyQ in  flowmatQ.keys():
            key="dilution"
            medium = fkeyQ[0]-1
            for entry in flowmat[fkeyQ]:
                comp = int(entry[0])-1
                for i in range (2,14):
                    season=i-2
                    if key in keylist[season][comp][medium].keys():
                        keylist[season][comp][medium][(key,entry[1])] += entry[i]
                    else:
                        keylist[season][comp][medium][(key,entry[1])] = entry[i]
    else:
        for fkey in  flowmat.keys():
            key="Flux"
            medium = fkey[0]-1
            for entry in flowmat[fkey]:
                comp = int(entry[0])-1
                # Remove all fluxes which are dealt with in Dproc            
                if entry[0] == entry[1]:
                    continue
                for i in range (2,14):
                    season=i-2
                    if key in keylist[season][comp][medium].keys():
                        keylist[season][comp][medium][(key,entry[1])] += entry[i]
                    else:
                        keylist[season][comp][medium][(key,entry[1])] = entry[i]
        

    return keylist


def mkDdiag_from_Fluxkey(fluxkey):
    Ddiag = sp.lil_matrix((12,lat*lon*10))
    for i1 in range(len(fluxkey)):
        for i2 in range(len(fluxkey[i1])):
            for i3 in range(len(fluxkey[i1][i2])):
                for key in fluxkey[i1][i2][i3].keys():
                    Ddiag[i1,i2*11+i3] += fluxkey[i1][i2][i3][key] 
    return Ddiag    



                    
def  match_keys(pkey):
        medium =  pkey[0]-1
        mstring = pkey[2]
        if re.search("deg", mstring):
            key = "deg"
        elif re.search("leach|burial|sink|strato|sedimentation", mstring):
            key = "sed"
        elif re.search("dep|diff|runoff|uptake|particle|disso|susp|mix|litter|erosion|sedi|exchange", mstring)  and (not mstring == "sedimentation"):
            key = pkey[1]-1
        else:
            raise KeyError(pkey[2])
        return key, medium

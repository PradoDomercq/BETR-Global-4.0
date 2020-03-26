################################################################################
'''Module "mkflowsD", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates D-values for inter-cell transport processes of a
particular model parametrization.'''
################################################################################
from numpy import *
import inspect
from globalz import *
import sys
import copy
import pdb

def mkflowD(model):
    fdict=copy.deepcopy(model.flowdict)
    #If velocity of wind G/A is smaller than 1000m/h, the G/A should be 1000m/h. here, we consider
    #upper and lower air.   by FY
    for [k,v] in fdict.items():
        if k[0]==1 and k[1]==1:
            for i in range (len(v)):
                for j in range(12):
                    a=v[i][0]-1
                    velocity=v[i][j+2]/(sqrt(model.par['A'][a][j])*model.par['h1'][a][j])
                    if velocity<1000:
                       v[i][j+2]=1000*(sqrt(model.par['A'][a][j])*model.par['h1'][a][j])
        if k[0]==2 and k[1]==2:
            for i in range (len(v)):
                for j in range(12):
                    a=v[i][0]-1
                    velocity=v[i][j+2]/(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])
                    if velocity<1000:
                       v[i][j+2]=1000*(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])

        if k[0]==8 and k[1]==8:
            for i in range (len(v)):
                for j in range(12):
                    a=v[i][0]-1
                    velocity=(v[i][j+2]*model.par['fcp1'][a][j])/(sqrt(model.par['A'][a][j])*model.par['h1'][a][j])
                    if velocity<1000:
                       v[i][j+2]=1000*(sqrt(model.par['A'][a][j])*model.par['h1'][a][j])
        if k[0]==9 and k[1]==9:
            for i in range (len(v)):
                for j in range(12):
                    a=v[i][0]-1 
                    velocity=(v[i][j+2]*model.par['ffp2'][a][j])/(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])                    
                    if velocity<1000:
                       v[i][j+2]=1000*(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])/model.par['ffp2'][a][j]       
        if k[0]==10 and k[1]==10:
            for i in range (len(v)):
                for j in range(12):
                    a=v[i][0]-1
                    velocity=(v[i][j+2]*model.par['fcp2'][a][j])/(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])
                    if velocity<1000:
                       v[i][j+2]=1000*(sqrt(model.par['A'][a][j])*model.par['h2'][a][j])/model.par['ffp2'][a][j]

    #Dflow=GZ
    for f in fdict.keys():        
        fromcells=fdict[f][:,0].astype(int)
        zvals=model.zdict[f[0]]['bulk'][fromcells-1,:]
        vvals=model.vdict[f[0]]['bulk'][fromcells-1,:]   #FY
        #fdict[f][:,2:]=fdict[f][:,2:]*zvals             #Old version
        # limited time is 2h   by FY
        time=6
        fdict[f][:,2:]=minimum(fdict[f][:,2:]*zvals, vvals*zvals/time)
    return(fdict)

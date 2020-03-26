################################################################################
'''Module "mkflowsDQ", October 2016, is part of
BETR-Research by Fangyuan Zhao

This module calculates D-values with Q for inter-cell transport processes of a
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

def mkflowDQ(m):       #Fy   15.08.2016
    comp=m.nocomp;     cells=m.nocells;  seasons=m.nots
    vz=zeros((comp,cells,seasons)); dk=zeros((comp,cells,seasons)) 
    Kloss=zeros((comp,cells,seasons))    #  Fy  20.09.2016
               
    ''' uses zdict and vdict to construct vz'''
    for c in m.compdict.keys():
        vz[c-1]=m.zdict[c]['bulk']*m.vdict[c]['bulk']
        
    '''use Dflow and Dprocess to construct Dtotalloss'''
    '''Dflow'''
    fdict=copy.deepcopy(m.Dflow)
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
        
    T=[]
    Q=[]
    for [k,v] in fdict.items():
        T.append([])
        Q.append([])
    for [k,v] in fdict.items():
        if k[0]==4 and k[1]==5: continue
        c=k[0]-1
        T[c]=zeros((len(v),len(v[0])))
        Q[c]=zeros((len(v),len(v[0])))
        for i in range(len(v)):
            a=v[i,0]-1
            T[c][i][0:2]=v[i,0:2]
            for ii in range(12):
                if v[i][(2+ii)]==0:
                   T[c][i][(2+ii)]=100000  #T can be defined as any values, the results won't be changed. because we only calculate upper and lower air
                else:
                   T[c][i][(2+ii)]=vz[c][a][ii]/v[i][(2+ii)]

    '''Dprocess'''  
    for [k,v] in m.Dproc.items():
        key=k[0]-1
        mstring = k[2]
        if re.search("deg|wet|disso", mstring): 
           if key==0: dk[0]+=v
           if key==1: dk[1]+=v
           if key==2: dk[2]+=v
           if key==3: dk[3]+=v
           if key==4: dk[4]+=v
           if key==5: dk[5]+=v
           if key==6: dk[6]+=v
           
    #here,we only redefine K in upper and lower air, and assume the K in the other compartments are 0.
    #so, the Q1,Q2 are calculated, and Q3-Q6 should be 1.
    for i in range(2):       
        for ii in range(cells):
            for iii in range(seasons):
                Kloss[i][ii][iii]=dk[i][ii][iii]/vz[i][ii][iii]
                
   
    '''calcuate Q, here we only calculate upper air and lower air, Q1,Q2'''
    for [k,v] in fdict.items():
        c=k[0]-1
        for ii in range(len(T[c])):
            b=T[c][ii][0]-1
            Q[c][ii][0:2]=T[c][ii][0:2]
            Q[c][ii][2:14]=(1+Kloss[c][b]*T[c][ii][2:14])/exp(Kloss[c][b]*T[c][ii][2:14])

    QQ={}
    QQ[0]=Q[0]; QQ[1]=Q[1]

    return(QQ)


def mkDQnocells(m):       #Fy   15.08.2016
    '''calculate total Q in every cell,
       such as for br, len(Q)=288'''
    
    comp=m.nocomp;     cells=m.nocells;  seasons=m.nots
    vz=zeros((comp,cells,seasons)); df=zeros((comp,cells,seasons)); dk=zeros((comp,cells,seasons)) 
    T=zeros((comp,cells,seasons)); K=zeros((comp,cells,seasons))    #  Fy  20.09.2016
    Q=zeros((comp,cells,seasons))
               
    ''' uses zdict and vdict to construct vz'''
    for c in m.compdict.keys():
        vz[c-1]=m.zdict[c]['bulk']*m.vdict[c]['bulk']
        
    '''use Dflow and Dprocess to construct Dtotalloss'''
    '''Dflow'''
    fdict=copy.deepcopy(m.Dflow)
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
        
    for [k,v] in fdict.items():
        mdf=zeros((comp,len(v[:,0]),14))
        if k[0]==4 and k[1]==5: continue
        c=k[0]-1
        mdf[c]=v
        for j in range(cells):
            for i in range(len(mdf[c])):
                if mdf[c][i][0]==(j+1):
                   df[c][j]+=mdf[c][i][2:14]   #c is compartment, j is fromcell,[2:14] are seasons.      

    for i in range(comp):
        for ii in range(cells):
            for iii in range(seasons):
                if df[i][ii][iii]==0:
                   T[i][ii][iii]=100000
                else:
                   T[i][ii][iii]=vz[i][ii][iii]/df[i][ii][iii]


    '''Dprocess'''
    for [k,v] in m.Dproc.items():
        key=k[0]-1
        if key==0: dk[0]+=v
        if key==1: dk[1]+=v
        if key==2: dk[2]+=v
        if key==3: dk[3]+=v
        if key==4: dk[4]+=v
        if key==5: dk[5]+=v
        if key==6: dk[6]+=v

    for i in range(comp):
        for ii in range(cells):
            for iii in range(seasons):
                if vz[i][ii][iii]==0:
                   K[i][ii][iii]=100000
                else:
                   K[i][ii][iii]=dk[i][ii][iii]/vz[i][ii][iii]
                   
   
    '''calcuate Q'''
    Q=(1+K*T)/exp(K*T)
    return (Q)
####################################################################################################

'''
    #Fangyuan 28th Nov, 2016
    #make Q in the sparse matrices form, and then we can mulitiply Q with Mass in solve.py
    Q=(1+K*T)/exp(K*T)
    # construct list of large sparse matrices
    Qlist=[]
    diags=zeros((m.nots,m.matdim))
    for i in range(m.nots):
        for ii in range(m.matdim):
            diags[i][ii]=1
    for c in range(2):
        idx=tocell(arange(1, m.nocells+1),(c+1),m)
        d=Q[c]
        for t in arange(0,m.nots):
            diags[t,idx]=d[:,t]
    for d in diags:
        v=compressvec(d,m.killidx)
        Qlist.append(sp.csr_matrix(v))
    return(Qlist)
'''

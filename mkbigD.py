################################################################################
'''Module "mkbigD.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

The function mkflowmat(m) takes a model object (m) as argument and returns
a list of sparse matrices which contain D-values from inter-cell flows
for each timestep.

The function mkprocmat(m) takes a model object (m) as argument and returns
a list of sparse matrices which contain D-values from intra-cell flows
for each timestep.

The function mkbigD(m) is called from BETRS.py, calls mkflowmat and mkprocmat
and returns the final list of system-matrices'''
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

def mkprocmat(m):
    ''' returns a list of sparse matrices, one for each timestep, that
    contain D-values for intra-cell processes'''  
    for tp in m.Dproc.keys(): # check whether we have all compartments
        tp=(tp[0],tp[1])
        for c in tp:
            if c not in m.compdict.keys():
                print('''mkprocmat.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # construct list of large sparse matrices
    matlist=[]
    #vz=zeros((m.nocomp,m.nocells,12))
    #for c in m.compdict.keys():
    #   vz[c-1]=m.zdict[c]['bulk']*m.vdict[c]['bulk']
        
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float32))))    #original: float64 
    for [k,v] in m.Dproc.items():
        toidx=array([x*m.nocomp+k[1]-1 for x in range(0,m.nocells)])
        fromidx=array([x*m.nocomp+k[0]-1 for x in range(0,m.nocells)])
        ij=(toidx.astype(int),fromidx.astype(int))
        
        count=0
        for ts in arange(0,v.shape[1]):
            matlist[count]=matlist[count]\
                         +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))
            count+=1
    return(matlist)

def mkprocmatQ(m):
    ''' returns a list of sparse matrices, one for each timestep, that
    contain D-values for intra-cell processes, consider mix process with Qmix,  Fy 14,06,2017'''
    
    pdict=copy.deepcopy(m.Dproc)  
    mixQ=copy.deepcopy(m.DmixQ)
    for tp in pdict.keys(): # check whether we have all compartments
        tp=(tp[0],tp[1])
        for c in tp:
            if c not in m.compdict.keys():
                print('''mkprocmat.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)

    for [k,v] in pdict.items():
        key=k[0]-1
        mstring = k[2]
        if re.search("mix", mstring): 
           if key==0:  v[:]=v[:]*mixQ[1][:]
           if key==1:  v[:]=v[:]*mixQ[0][:]

    # construct list of large sparse matrices
    matlist=[]
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float32))))   #original: float64 
    for [k,v] in pdict.items():
        toidx=array([x*m.nocomp+k[1]-1 for x in range(0,m.nocells)])
        fromidx=array([x*m.nocomp+k[0]-1 for x in range(0,m.nocells)])
        ij=(toidx.astype(int),fromidx.astype(int))
        count=0
        for ts in arange(0,v.shape[1]):
            matlist[count]=matlist[count]\
                         +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))
            count+=1          
    
    return(matlist)

def mkflowQ(m):
    '''calcuate Gn-1Qn
    Correct the inflow rate with Q '''
    
    fdict=copy.deepcopy(m.Dflow)
    Qdict=copy.deepcopy(m.DflowQ)      #04.10.2016   fy
    for tp in fdict.keys(): # check whether we have all compartments
        for c in tp:
            if c not in m.compdict.keys():
                print('''flows.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # remove intra-cell flows
    # this is for example oceanic sinking flux, which is dealt with
    # by "betr_ocean_sinkflux" in "processes.py".
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
     # construct list of large sparse matrices

    #GN-1Qn
    col1=[];   col2=[]
    if m.nocells==288:
       raw0=12;  col0=24;  da1=24;   da2=23       
    if m.nocells==1152:
       raw0=24;  col0=48;  da1=48;   da2=47 
    if m.nocells==4608:
       raw0=48;  col0=96;  da1=96;   da2=95  

    for i in arange(0,raw0):
        x=col0*i+1
        y=col0*(i+1)
        col1.append(x)
        col2.append(y)
    
    for [kk,vv] in fdict.items():
        if len(vv)==0: continue
        for [k,v] in Qdict.items():
            if (kk[0]-1)==k:
               nv=zeros((len(v),14))           
               for ii in range(len(v)):
                   diffv=v[ii,1]-v[ii,0]                
                   nv[ii,0]=v[ii,0];  nv[ii,1]=v[ii,1]
                   if diffv==da1 or diffv==-da1:        
                       a=v[ii,1]; b=v[ii,1]+diffv
                       for iii in range(len(v)):
                           if v[iii,0]==a and v[iii,1]==b:
                              nv[ii,2:14]=v[iii,2:14]
                   if diffv==1:
                      a=v[ii,1]
                      if a in col2:
                         b=v[ii,1]-da2     
                      else:
                         b=v[ii,1]+diffv
                      for iii in range(len(v)):
                          if v[iii,0]==a and v[iii,1]==b:
                             nv[ii,2:14]=v[iii,2:14]              
                   if diffv==-1:
                      a=v[ii,1]
                      if a in col1:
                         b=v[ii,1]+da2     
                      else:
                         b=v[ii,1]+diffv
                      for iii in range(len(v)):
                          if v[iii,0]==a and v[iii,1]==b:
                             nv[ii,2:14]=v[iii,2:14]
                   if diffv==da2:          
                      a=v[ii,1]; b=v[ii,1]-1
                      for iii in range(len(v)):
                         if v[iii,0]==a and v[iii,1]==b:
                            nv[ii,2:14]=v[iii,2:14]
                   if diffv==-da2:     
                      a=v[ii,1]; b=v[ii,1]+1
                      for iii in range(len(v)):
                          if v[iii,0]==a and v[iii,1]==b:
                              nv[ii,2:14]=v[iii,2:14]

               for i in range(len(v)):
                    if nv[i,2]==0 and nv[i,3]==0 and nv[i,4]==0 and nv[i,5]==0 and nv[i,6]==0 and nv[i,7]==0\
                       and nv[i,8]==0 and nv[i,9]==0 and nv[i,10]==0 and nv[i,11]==0\
                       and nv[i,12]==0 and nv[i,13]==0:
                       anocell=nv[i,1]-1;   nok=k
                       nv[i,2:14]=calQ(nok,anocell,m)
                       #nv[i]=v[i]
               vv[:,2:14]=vv[:,2:14]*v[:,2:14]
    return(fdict)

def mkflowmat(m):
    '''takes a model object (m) as argument and returns a list of sparse
    matrices which contain D-values from inter-cell flows for each timestep.'''
    fdict=copy.deepcopy(m.Dflow)
    for tp in fdict.keys(): # check whether we have all compartments
        for c in tp:
            if c not in m.compdict.keys():
                print('''flows.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # remove intra-cell flows
    # this is for example oceanic sinking flux, which is dealt with
    # by "betr_ocean_sinkflux" in "processes.py".
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
     # construct list of large sparse matrices
    matlist=[]
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float32))))  # original: float64 
    for [k,v] in fdict.items():
        if len(v)==0: continue
        toidx=array([(x-1)*m.nocomp+k[1]-1 for x in v[:,1]])
        fromidx=array([(x-1)*m.nocomp+k[0]-1 for x in v[:,0]])
        ij=(toidx.astype(int),fromidx.astype(int))
        count=0
        for ts in arange(2,v.shape[1]):
            matlist[count]=matlist[count]\
                        +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))           
            count+=1
    return(matlist)

def mkflowmatQ(m):
    '''Return a list of sparse
    matrices which contain D-values with Q from inter-cell flows for each timestep.'''
    fdict=mkflowQ(m)
    for tp in fdict.keys(): # check whether we have all compartments
        for c in tp:
            if c not in m.compdict.keys():
                print('''flows.py: Compartment %i not in compartment-list.
                Aborting !''') % (c)
                sys.exit(1)
    # remove intra-cell flows
    # this is for example oceanic sinking flux, which is dealt with
    # by "betr_ocean_sinkflux" in "processes.py".
    for [k,v] in fdict.items():
        intercell=list(where(v[:,0]!=v[:,1])[0])
        intracell=list(where(v[:,0]==v[:,1])[0])
        fdict[k]=v[intercell,:]
     # construct list of large sparse matrices

    matlist=[]
    for t in arange(0,m.nots):
        matlist.append(sp.coo_matrix(matrix(zeros((m.matdim,m.matdim),
                                                  dtype=float32))))  # original: float64 
    for [k,v] in fdict.items():
        if len(v)==0: continue
        toidx=array([(x-1)*m.nocomp+k[1]-1 for x in v[:,1]])
        fromidx=array([(x-1)*m.nocomp+k[0]-1 for x in v[:,0]])
        ij=(toidx.astype(int),fromidx.astype(int))
        count=0
        for ts in arange(2,v.shape[1]):
            matlist[count]=matlist[count]\
                        +sp.coo_matrix((v[:,ts],ij),shape=(m.matdim,m.matdim))
            count+=1
    return(matlist)

def mkbigD(m, track_flows=False ,use_correction=False):
    '''with Q correction'''    
    matlist_flow=mkflowmat(m)
    matlist_proc=mkprocmat(m)
    if use_correction==True:
        matlist_flow_Q=mkflowmatQ(m)
        matlist_proc_Q=mkprocmatQ(m)

    # check for same shape of all matrices
    matshp=matlist_flow[0].shape
    for m in matlist_flow+matlist_proc:
        if m.shape != matshp:
            print("mkbigD: inconsistent matrix dimensions. Aborting!\n")
            sys.exit(1)
	
    # add matrices for intercell transport  and intracell processes 
    matlist=[]; outbigD=[]
    if use_correction==True:
       matlist_Q=[]
    for ts in arange(0,len(matlist_flow)):
        matlist.append(matlist_flow[ts]+matlist_proc[ts])
        outbigD.append(matlist_flow[ts]+matlist_proc[ts])
        if use_correction==True:
           matlist_Q.append(matlist_flow_Q[ts]+matlist_proc_Q[ts])
	
    # construct proper diagonals
    for ts in arange(0,len(matlist)):
        # seperate matrices from their diagonals        
        # define todense() since it can not work in high resolution   FY
        array_matlist1=zeros((matlist[ts].shape[0],matlist[ts].shape[1]))       
        for i in range(matlist[ts].shape[0]):
            array_matlist1[i,:]=matlist[ts][i,:].todense()
        #diagonal=diag(matlist[ts].todense())
        diagonal=diag(array_matlist1)        #FY
        matlist[ts]=matlist[ts]-sp.spdiags(diagonal,0,matshp[0],matshp[1],
                                        format="csc")
        if use_correction==True:
            diagonal_Q=diag(matlist_Q[ts].todense())      
            matlist_Q[ts]=matlist_Q[ts]-sp.spdiags(diagonal_Q,0,matshp[0],matshp[1],
                                        format="csc")
        # put losses by transport to other cells on diagonal and
        # put back intra-cell loss with negative sign
        #loss=sum(matlist[ts].todense(), axis=0)
        array_matlist2=zeros((matlist[ts].shape[0],matlist[ts].shape[1]))       
        for i in range(matlist[ts].shape[0]):
            array_matlist2[i,:]=matlist[ts][i,:].todense()
        loss=sum(array_matlist2, axis=0)   #FY
        if use_correction==True:
            outbigD[ts]=matlist_Q[ts]\
                       -sp.spdiags(loss,0,matshp[0],matshp[1],format="csc")\
                       -sp.spdiags(diagonal,0,matshp[0],matshp[1],format="csc")
        else:
            outbigD[ts]=matlist[ts]\
                       -sp.spdiags(loss,0,matshp[0],matshp[1],format="csc")\
                       -sp.spdiags(diagonal,0,matshp[0],matshp[1],format="csc")                     
        if track_flows==True:
            diagonals = -sp.spdiags(loss,0,matshp[0],matshp[1],format="csc")\
                    -sp.spdiags(diagonal,0,matshp[0],matshp[1],format="csc")
            outbigD[ts] = sp.bmat([[ outbigD[ts], None],[diagonals,\
            sp.csc_matrix(numpy.zeros(matshp))]])
    return(outbigD)

def mkZVinv(m):
    ''' uses zdict and vdict to construct a list of
    sparse (csr) diagonal matrices A = 1/ZV'''
    # construct list of large sparse matrices
    mlist=[]
    diags=zeros((m.nots,m.matdim))
    for c in m.compdict.keys():
        idx=tocell(arange(1, m.nocells+1),c,m)
        d=m.zdict[c]['bulk']*m.vdict[c]['bulk']    
        for t in arange(0,m.nots):
            diags[t,idx]=d[:,t]
    for d in diags:
        v=1/compressvec(d,m.killidx)
        mlist.append(sp.spdiags(v,0,len(v),len(v), format='csr'))
    return(mlist)

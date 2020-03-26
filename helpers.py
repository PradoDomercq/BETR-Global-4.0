################################################################################
'''Module "helpers.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module contains various helper functions'''
################################################################################
from numpy import *
from globalz import *
import scipy.sparse as sp
import sys
import re
import pdb

def findemptycells_original(A): # This is the original findemptycells. RKG. 10.06.2014
    """ returns a list of indices with zero rows resp. columns"""
    B=A.toarray()
    idx=[]
    for i in range(0,B.shape[0]):
        if (all(B[i,:] == 0)):
            idx.append(i)
            if (not all(B[:,i] == 0)):
                sys.exit("findemptycells: empty cell receiving some ! Abort !")
    return(idx)

def findemptycells(A):  # function added by RKG, 10.06.2014
                            # This is an alternative to findemptycells (see above) that does the same thing...
                            # ...without using the toarray() method on the whole matrix...
                            # ...but on individual rows and columns. RKG
    """ returns a list of indices with zero rows resp. columns"""
    #print('findemptycells: returns a list of indices with zero rows resp. columns')
    idx=[]
    for i in range(0,A.shape[0]):
        if(all(A.getrow(i).toarray() == 0)):
            idx.append(i)
            if (not all(A.getcol(i).toarray() == 0)):
                sys.exit("findemptycells: empty cell receiving some ! Abort !")
    return(idx)

def findnonemptycells(A):      # function added by HW, 17.01.2011
    """ returns a list of indices with zero rows resp. columns"""
    idx=[]
    for i in range(0,A.shape[0]):
        if (any(A[i,:] != 0)):
            idx.append(i)
    return(idx)

def compressmat(A, killidx):
    ''' deletes rows and columns in killidx from matrix'''
    spformat=A.getformat()
    dim= shape(A)[0]
    keepidx =[i for i in range(shape(A)[0]) if not i in killidx]
    X = sp.coo_matrix((ones(len(keepidx)),(keepidx, arange(0,len(keepidx)))), shape = (dim, len(keepidx)))
    return (X.T * A * X).asformat(spformat)

def compressvec(v, killidx):
    """ Takes vector v and returns weeded vector"""
    vw = delete(v, killidx)
    return vw

def expandmatrix(A, killidx):
    spformat=A.getformat()
    '''re-expands previosly compressed matrix'''
    killidx=sort(killidx)
    B=A.toarray()
    for i in killidx:
        B = insert(B,i,[0],axis=0)
        B = insert(B,i,[0],axis=1)
    B = sp.coo_matrix(matrix(B))
    return(B.asformat(spformat))

def expandvec(v, killidx):
    """ Re - inserts zeros at the killidx - positions"""
    ve=copy(v)
    killidx=sort(killidx)
    for i in killidx:
        ve=insert(ve,i,0)
    return ve

def cell2regcomp(cell,m):
    '''  translates cell index into region-number and compartment ID '''
    reg=floor(cell/m.nocomp) + 1
    comp=mod(cell,m.nocomp) + 1
    return (reg.astype(int),comp.astype(int))

def tocell(reg,comp,m):
    ''' translates region number and compartment ID into cell-index'''
    cell=m.nocomp*(reg-1)+comp-1 
    return cell

def calQ(comp,cell,m):
    dk=zeros((2,m.nocells,12))
    x=zeros(12)
    for [k,v] in m.Dproc.items():
        key=k[0]-1
        mstring = k[2]
        if re.search("deg|wet|disso", mstring): 
           if key==0: dk[0]+=v
           if key==1: dk[1]+=v
           
    if comp==0:
        for j in range(12):
            x[j]=1000*(sqrt(m.par['A'][cell][j])*m.par['h1'][cell][j])
        KT=dk[0][cell]/(m.zdict[1]['bulk'][cell]*x)
        Q=(1+KT)/exp(KT)
        
    if comp==1:
        for j in range(12):
            x[j]=1000*(sqrt(m.par['A'][cell][j])*m.par['h2'][cell][j])
        KT=dk[1][cell]/(m.zdict[2]['bulk'][cell]*x)
        Q=(1+KT)/exp(KT)
    return(Q)
  
def calGZ(comp,cell,m):
    x=zeros(12)           
    if comp==0:
        for j in range(12):
            x[j]=1000*(sqrt(m.par['A'][cell][j])*m.par['h1'][cell][j])
        GZ=m.zdict[1]['bulk'][cell]*x
    if comp==1:
        for j in range(12):
            x[j]=1000*(sqrt(m.par['A'][cell][j])*m.par['h2'][cell][j])
        GZ=m.zdict[2]['bulk'][cell]*x

    return(GZ)
     

###############################################################################

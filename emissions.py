################################################################################
'''Module "emissions.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module reads emission files and hands out emission-vectors'''
################################################################################
from numpy import *
import scipy.sparse as sp
import inspect
from globalz import *
from helpers import *
import sys
import pdb
    
class Emission:
    """ reads emission file and saves as static dictionary
        The compartments start with 1 """  
       
    def __init__(self,fn):
        self.em = {}
        try:
            f=open(fn, 'r')
        except IOError:
            print("emissions.py: emission file %s not found. Aborting!") % (fn)
            sys.exit(1)
        lines=f.readlines()
        f.close
        for lin in lines:
            if (lin[0] == '#' or lin == ''):
                continue
            #Fangyuan changed [m,r,c,val] to [m,r,c1,c2,c3,val,perc1,perc2]
            #Fangyuan changed append((int(float(r)),int(float(c)),float(val))) to
            #append((int(float(r)),int(float(c1)),int(float(c2)),int(float(c3)),float(val)),int(float(perc1)),int(float(perc2)))
            [m,r,c1,c2,c3,val,perc1,perc2]=[x.lstrip() for x in \
                         [x.rstrip() for x in lin.split()]]
            self.em.setdefault(int(float(m)),[]).\
                             append((int(float(r)),int(float(c1)),int(float(c2)),int(float(c3)),
                                     float(val),float(perc1),float(perc2)))
            
    def get_emission(self,timeidx,m,type='csc'):
        """ hands out emission vector for timeindex (starts with 0).
        default type csc returns vector in csc-format
        type=array returns 1D-array."""
        emvec=zeros(m.matdim)
        try:
            ems=self.em[timeidx+1]
            for source in ems:
                #Fangyuan separate emissions to gas phase, fine and coarse particles. 
                cell1 = tocell(source[0],source[1],m)
                emvec[cell1]=source[4]*source[5]                
                cell2 = tocell(source[0],source[2],m)
                emvec[cell2]=source[4]*(1-source[5])*source[6]                
                cell3 = tocell(source[0],source[3],m)
                emvec[cell3]=source[4]*(1-source[5])*(1-source[6])                
        except KeyError:
            pass
        ## compress
        emvec=compressvec(emvec,m.killidx)
        if type=='csc':
             emvec=sp.csr_matrix(emvec).T
        return emvec
    
    # def get_csc(self,timeidx,m):
    #     """ hands out emission vector for timeindex (starts with 0) as csc-matrix"""
    #     emvec=zeros(m.matdim)
    #     try:
    #         ems=self.em[timeidx+1]
    #         for source in ems:
    #             cell = tocell(source[0],source[1],m)
    #             emvec[cell]=source[2]
    #     except KeyError:
    #         pass
    #     ## compress and change to csc-vector
    #     emvec=sp.csr_matrix(compressvec(emvec,m.killidx)).T
    #     # emvec=compressvec(emvec,m.killidx)
    #     ## change units of mol/h to fugaity [Pa/h]
    #     # season=mod(timeidx,m.nots)
    #     # emvec=(m.ZVinvlist[season]*sp.csr_matrix(emvec).T).asformat('csc')
    #     return emvec


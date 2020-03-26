################################################################################
'''Module "volumes", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates compartments volumes'''
################################################################################
from numpy import *
import inspect
from globalz import *

class Volumes():
    '''This class constructs a dictionary with compartment-IDs as keys and
    dictionaries with bulk- and sub-compartment volumes as values.
    It calls all a method V<compartment_ID> for each compartment in compdict.
    This method has to exist, and must return a dictionary with keys "bulk" and
    optionally sub-compartment names.

    In case no volumes are necessary to calculate processes for a particular
    compartment with ID <x>, or volumes are directly read from an input-file,
    define a function V<x> that returns an empty dictionary.'''

## dont't touch this ###########################################################
    def __init__(self, par, compdict):
        self.par=par
        self.compdict=compdict
        empty=[]
        for i in arange(0,len(compdict)):
            empty.append({})
        self.V=dict(zip(self.compdict.keys(), empty))
        for c in self.compdict.keys():
            try:
                self.V[c]=getattr(self,'V'+str(c))()
            except AttributeError:
                sys.exit('Calculation of volumes for compartment {0:d} '
                         +'not implemented. Aborting !'.format(c))
            
################################################################################
    
## DEFINE VOLUMES FOR NEW COMPARTMENTS HERE ####################################
    def V1(self):   #FY
        vdict={}
        vdict['bulk']=self.par['h1']*self.par['A']*(1-self.par['fcp1'])
        return(vdict)
    
    def V2(self):     #FY
        vdict={}
        vdict['bulk']=self.par['h2']*self.par['A']*(1-self.par['ffp2']-self.par['fcp2'])
        return(vdict)
    
    def V3(self):
        vdict={}
        rhoveg=self.par['fw3']*self.par['rho45']\
               +(1-self.par['fw3'])*self.par['rho3']
        vdict['bulk']=self.par['A']*self.par['perc6']*self.par['perc3']\
                      *self.par['mveg']/rhoveg 
        vdict['water']=vdict['bulk']*self.par['fw3']
        vdict['flesh']=vdict['bulk']*(1-self.par['fw3'])
        return(vdict)

    def V4(self):
        vdict={}
        vdict['bulk']=self.par['A']*self.par['perc4']*self.par['h4']
        vdict['sussed']=vdict['bulk']*self.par['fp4']
        vdict['biota']=vdict['bulk']*self.par['ff4']
        vdict['water']=vdict['bulk']-vdict['sussed']-vdict['biota']
        return(vdict)
    
    def V5(self):
        vdict={}
        vdict['bulk']=self.par['A']*self.par['perc5']*self.par['h5']
        vdict['sussed']=vdict['bulk']*self.par['fp5']
        vdict['biota']=vdict['bulk']*self.par['ff5']
        vdict['water']=vdict['bulk']-vdict['sussed']-vdict['biota']
        return(vdict)
    
    def V6(self):
        vdict={}
        vdict['bulk']=self.par['A']*self.par['perc6']\
                      *self.par['h6']
        vdict['air']=vdict['bulk']*self.par['fa6']
        vdict['water']=vdict['bulk']*self.par['fw6']
        vdict['solids']=vdict['bulk']*self.par['fs6']
        return(vdict)

    def V7(self):
        vdict={}
        vdict['bulk']=self.par['A']*self.par['perc4']*self.par['h7']
        vdict['water']=vdict['bulk']*self.par['fw7']
        vdict['solids']=vdict['bulk']*self.par['fs7']
        return(vdict)

    def V8(self):     #FY
        vdict={}
        vdict['bulk']=self.par['h1']*self.par['A']*self.par['fcp1']
        return(vdict)

    def V9(self):    #FY
        vdict={}
        vdict['bulk']=self.par['h2']*self.par['A']*self.par['ffp2']
        return(vdict)

    def V10(self):    #FY
        vdict={}
        vdict['bulk']=self.par['h2']*self.par['A']*self.par['fcp2']
        return(vdict)

###############################################################################


    
    
   
       




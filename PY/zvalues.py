################################################################################
'''Module "zvalues", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates Z-values'''
################################################################################

from numpy import *
import inspect
from globalz import *
import copy
import pdb

class Zvalues():
    '''This class constructs a dictionary with compartment-IDs as keys and
    dictionaries with bulk- and sub-compartment Z-values (mol/Pa/m^3) as values.
    It calls all a method Z<compartment_ID> for each compartment in compdict.
    This method has to exist, and must return a dictionary with keys "bulk" and
    optionally sub-compartment names.

    In case no Z-values are necessary to calculate processes for a particular
    compartment with ID <x>, define a function Z<x> that returns an empty
    dictionary.'''

## dont't touch this ###########################################################
    def __init__(self, par, compdict, chempardict):
        self.par=par
        self.compdict=compdict
        self.chempardict=chempardict
        empty=[]
        for i in arange(0,len(compdict)):
            empty.append({})
        self.Z=dict(zip(self.compdict.keys(), empty))       
        for c in self.compdict.keys():
            try:
                self.Z[c]=getattr(self,'Z'+str(c))()
            except AttributeError:
                sys.exit('Calculation of Z-values for compartment {0:d} '
                         +'not implemented. Aborting !'.format(c))
            
################################################################################
    
## DEFINE Z-VALUES FOR NEW COMPARTMENTS HERE ###################################
#
#    Functions that calculate Z-values for different compartments
#    and sub-compartments. Each function returns a dictionary containing
#    the sub-compartment names and/or 'bulk' as keys, and the
#    Z-values (mol/Pa/m^3) as values.
    
    def Z1(self):
        '''Z-values for upper air'''
        T=self.par[self.compdict[1]['temp_variable']] #temperature field
        zdict={}
        zdict['bulk']=(R*T)**-1        
        return(zdict)

    def Z2(self):
        '''Z-values for lower air'''
        T=self.par[self.compdict[2]['temp_variable']] #temperature field
        zdict={}
        zdict['bulk']=1/(R*T)
        ## this is a pseudo Z-value that takes into account the
        ## changed bulk Z-value of the atmosphere during rain events.
        ### ## volume fraction of water in air during rain
        norainmask = self.par['stwet'] == 0      
        stwet_tmp=copy.copy(self.par['stwet'])
        stwet_tmp[norainmask]=1
        mtc_event=self.par['precip']\
                   *(self.par['stdry']+self.par['stwet'])/stwet_tmp
        mtc_event[norainmask]=0
        fvw=mtc_event/self.par['vrain']
        zrain=zdict['bulk']/self.chempardict[2]['Kaw'] # actual Z-value of rain
        zdict['rain']=zrain*zdict['bulk']/(fvw*zrain+(1-fvw)*zdict['bulk'])   #DOUBLE CHECK!
        return(zdict)

    def Z3(self):
        '''Z-values for vegetation'''
        T=self.par[self.compdict[3]['temp_variable']] #temperature field
        zdict={}
        zdict['water']=(R*T*self.chempardict[3]['Kaw'])**-1
        # vegetation-flesh - water partitioning coefficient.
        # fo3 = volume fraction of pseudo octanol in vegetation-flesh
        Kvegw = self.par['fo3']*self.chempardict[3]['Kow']
        zdict['flesh']=zdict['water']*Kvegw
        zdict['bulk']=self.par['fw3']*zdict['water']+\
                       (1-self.par['fw3'])*zdict['flesh']
        return(zdict)

    def Z4(self):
        '''Z-values for freshwater'''
        T=self.par[self.compdict[4]['temp_variable']] #temperature field
        zdict={}
        zdict['water']=(R*T*self.chempardict[4]['Kaw'])**-1
        # Kqw: suspended sediment - water partitioning coeffcient
        Kqw=0.41*self.chempardict[4]['Kow'] # organic carbon - water [L/kg]
        Kqw*=self.par['focp4'] # susp. sediment - water [L/kg]      
        Kqw*=self.par['rhop45']/1000 #  susp. sediment - water [L/L]
        zdict['sussed']=zdict['water']*Kqw
        # Kfw fish(biota) - water partitioning coefficient
        # ATT : reference for '0.05' in Kfw=0.05*Kow ?
        Kfw=0.05*self.chempardict[4]['Kow']
        zdict['biota']=zdict['water']*Kfw
        zdict['bulk']=(1-self.par['fp4']-self.par['ff4'])*zdict['water']+\
                      self.par['fp4']*zdict['sussed']+\
                      self.par['ff4']*zdict['biota']
        # Add Z-value of air-boundary layer over fresh-water, with T = T_freshwater # added by HW
        zdict['air']=1/(R*T)
        return(zdict)

        
    def Z5(self):
        '''Z-values for ocean water'''
        T=self.par[self.compdict[5]['temp_variable']] #temperature field
        zdict={}
        # ATT : factor 0.8 in Zoceanwater=0.8*Zfreshwater from BETR-VBA. Why ?
        zdict['water']=0.8*(R*T*self.chempardict[5]['Kaw'])**-1
        # Kqw: suspended sediment - water partitioning coeffcient
        # ATT : replace 0.41 (from Karikoff 1981) by 0.35 (from Seth 1999)
        Kqw=0.41*self.chempardict[5]['Kow'] # organic carbon - water [L/kg]
        Kqw*=self.par['focp5'] # susp. sediment - water [L/kg]
        Kqw*=self.par['rhop45']/1000 #  susp.sediment - water [L/L]
        zdict['sussed']=zdict['water']*Kqw
        # Kfw fish(biota) - water partitioning coefficient
        # ATT : reference for '0.05' in Kfw=0.05*Kow ?
        Kfw=0.05*self.chempardict[5]['Kow']
        zdict['biota']=zdict['water']*Kfw
        zdict['bulk']=(1-self.par['fp5']-self.par['ff5'])*zdict['water']+\
                      self.par['fp5']*zdict['sussed']+\
                      self.par['ff5']*zdict['biota']
        # Add Z-value of air-boundary layer over ocean-water, with T = T_oceanwater # added by HW
        zdict['air']=1/(R*T)
        return(zdict)
        

    def Z6(self):
        '''Z-values for soil'''
        T=self.par[self.compdict[6]['temp_variable']] #temperature field
        zdict={}
        zdict['air']=1/(R*T)
        zdict['water']=zdict['air']/self.chempardict[6]['Kaw']
        # Ksw : solids in soil - water partitioning coefficient
        Ksw=0.41*self.chempardict[6]['Kow'] # organic carbon - water [L/kg]
        Ksw*=self.par['focs6'] # solid fraction - water [L/kg]
        Ksw*=self.par['rhos6']/1000 #  solid fraction - water [L/L]
        zdict['solids']=zdict['water']*Ksw
        zdict['bulk']=self.par['fa6']*zdict['air']+\
                      self.par['fw6']*zdict['water']+\
                      self.par['fs6']*zdict['solids']
        return(zdict)

    def Z7(self):
        '''Z-values for sediment'''
        T=self.par[self.compdict[7]['temp_variable']] #temperature field
        zdict={}
        zdict['water']=(R*T*self.chempardict[7]['Kaw'])**-1
        # Ksedw : solids in sediment - water partitioning coefficient
        Ksedw=0.41*self.chempardict[7]['Kow'] # organic carbon - water [L/kg]
        Ksedw*=self.par['focs7'] # solid fraction - water [L/kg]
        Ksedw*=self.par['rhos7']/1000 #  solid fraction - water [L/L]
        zdict['solids']=zdict['water']*Ksedw
        zdict['bulk']=self.par['fw7']*zdict['water']+\
                      self.par['fs7']*zdict['solids']
        return(zdict)

    def Z8(self):
        '''Z-values for coarse aerosol in upper air'''
        T=self.par[self.compdict[1]['temp_variable']] #temperature field
        zdict={}
        # aerosol-air partitioning coefficient  Kp=Fom*1500/840*0.26*Koa
        Kqa=self.par['cfom1']*1500/840*0.6*self.chempardict[1]['Koa']
        zdict['bulk']=(R*T)**-1*Kqa   
        return(zdict)

    def Z9(self):
        '''Z-values for fine aerosol in lower air'''
        T=self.par[self.compdict[2]['temp_variable']] #temperature field
        zdict={}
        # aerosol-air partitioning coefficient  Kp=Fom*1500/840*0.26*Koa
        Kqa=self.par['ffom2']*1500/840*0.6*self.chempardict[2]['Koa']
        zdict['bulk']=(R*T)**-1*Kqa   
        return(zdict)

    def Z10(self):
        '''Z-values for coarse aerosol in lower air'''
        T=self.par[self.compdict[2]['temp_variable']] #temperature field
        zdict={}
        # aerosol-air partitioning coefficient  Kp=Fom*1500/840*0.26*Koa
        Kqa=self.par['cfom2']*1500/840*0.6*self.chempardict[2]['Koa']
        zdict['bulk']=(R*T)**-1*Kqa   
        return(zdict)

###############################################################################


# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:23:01 2013
secondary_emission_reconstruction
@author: sebastian
"""
#from numpy import *
#from globalz import *

import scipy.sparse as sp
#import numpy as np
from numpy import zeros, shape, array, asarray

from helpers3 import mk_compmap, x_to_cellcomp, normalizefluxkey


def init_szenarios():
    pass

def reconstruct_dyn_res(m, scenario="air",nseasons=12):
    ''' Reconstruct the time resolved solution for masses from fluxes '''
    
    scenarios = {"air" : [[ 1, 1, 1, 1, 1, 1, 1], #air to *
                          [1, 1, 1, 1, 1, 1, 1],  #air to *
                          [0, 0, 1, 0, 0, 0, 0],  #veg to *
                          [0, 0, 0, 1, 1, 0, 0],  #wat to *
                          [0, 0, 0, 1, 1, 0, 0],  #wat to *
                          [0, 0, 0, 0, 0, 1, 0],  #soil to *
                          [0, 0, 0, 0, 0, 0, 1]], #sed to *
                          
                "water" : [[ 1, 1, 1, 0, 0, 1, 1],
                          [1, 1, 1, 0, 0, 1, 1],
                          [0, 0, 1, 0, 0, 0, 0],
                          [1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 0, 0, 1, 0],
                          [0, 0, 0, 0, 0, 0, 1]],
    
                "soil" : [[ 1, 1, 1, 1, 1, 0, 1],
                          [1, 1, 1, 1, 1, 0, 1],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [1, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0]],
    
                "sed" : [[ 1, 1, 1, 1, 1, 1, 0],
                          [1, 1, 1, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 1, 1, 1, 1]],
    
                "veg" : [[ 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0]],
                }
    
    inherits_mat = array(scenarios[scenario])
    
    
    dyn_res_rec = zeros(shape(m.dyn_res))
    dyn_res_pe = zeros(shape(m.dyn_res))
    dyn_res_se = zeros(shape(m.dyn_res))
    
    dyn_res_rec[:,0] = m.dyn_res[:,0]
    dyn_res_pe[:,0] = m.dyn_res[:,0]
    compmap = mk_compmap(m)     
    
    flux_key = normalizefluxkey(m.flux_key)
    
    for ts in range(len(m.flux_res)):
#        print "ts = %i" % (ts)
        '''Load massed from t-1 '''
        org_mass = m.dyn_res[:,ts]
        org_mass_pe = dyn_res_pe[:,ts]
        org_mass_se = dyn_res_se[:,ts]
        
        org_frac_se = calc_org_frac_se(org_mass_pe, org_mass_se)
                
        plus_mass_se = zeros(shape(org_mass))
        plus_mass_pe = zeros(shape(org_mass))
        minus_mass = zeros(shape(org_mass))
        
        season = ts % nseasons
        year = ts / nseasons
        
        fluxmat = m.flux_res[ts]
        x,y = fluxmat.nonzero()
        
        if season ==0:
            print ("\nProcessing year %i \t" %(year)),
        
        print ".",
        
        ''' For loop through compartment and treat in fluxes first '''
        for i in range(len(x)):
#            x,y = fluxmat.nonzero()
            value = fluxmat[x[i],y[i]]            
            xcell, xcomp = x_to_cellcomp( x[i], compmap, m.nocomp)
#            ycell, ycomp = x_to_cellcomp( y[i], compmap, m.nocomp)
            
            ind_from = xcell*m.nocomp + xcomp
#            ind_to = xcell*m.nocells + xcomp
            
            if x[i]==y[i]:
                assert value <= 0.
                try :
                    match_dict = m.flux_key[season][xcell][xcomp]
                except KeyError:
                    print [[season],[xcell],[xcomp]]
                
                minus_mass[ind_from] -= value
                
                for mnkey in match_dict.keys():
                     
                    '''Flow to and from neighboring cells'''
                    if type(mnkey) is tuple:
#                        minus_mass[ind_from] -= \
#                            flux_key[season][xcell][xcomp][mnkey] * \
#                            value
                        try: 
                            plus_mass_pe[(mnkey[1]-1)*m.nocomp+xcomp] -= \
                                flux_key[season][xcell][xcomp][mnkey] * \
                                value * \
                                (1-org_frac_se[ind_from])
                                
                            plus_mass_se[(mnkey[1]-1)*m.nocomp+xcomp] -= \
                                flux_key[season][xcell][xcomp][mnkey] * \
                                value * \
                                (org_frac_se[ind_from])
                        except IndexError:
                            print mnkey[1], (mnkey[1]-1)*m.nocomp+xcomp, ind_from
                            raise IndexError
#                    if mnkey  == "deg":
#                        minus_mass[ind_from] -= flux_key[season][xcell][xcomp][mnkey] * value
#                    if mnkey == "sed":
#                        minus_mass[ind_from] -= flux_key[season][xcell][xcomp][mnkey] * value
                    
                    if type(mnkey) == int:
#                        minus_mass[ind_from] -= \
#                            flux_key[season][xcell][xcomp][mnkey] * \
#                            value
                        
                        if inherits_mat[xcomp, mnkey] == 1:
                            
                            plus_mass_pe[xcell*m.nocomp+mnkey] -= \
                                flux_key[season][xcell][xcomp][mnkey] * \
                                value * \
                                (1-org_frac_se[ind_from])
                            
                            plus_mass_se[xcell*m.nocomp+mnkey] -= \
                                flux_key[season][xcell][xcomp][mnkey] * \
                                value * \
                                (org_frac_se[ind_from])
                                
                        if inherits_mat[xcomp, mnkey] == 0:
#                            print "HERE"
                            plus_mass_se[xcell*m.nocomp+mnkey] -= \
                                flux_key[season][xcell][xcomp][mnkey] * \
                                value
#                            print max(plus_mass_se)
        
        '''Calculate reconstructed mass from fluxes'''
        rec_mass =  org_mass + plus_mass_se + plus_mass_pe - minus_mass
        '''Calculate emissions from missing mass'''
        
        emissions = m.dyn_res[:,ts+1] - rec_mass
        
#        return (m.dyn_res[:,ts+1], emissions, org_mass, plus_mass_se, plus_mass_pe, minus_mass)
        '''Add emissions to primary emissions'''
        plus_mass_pe += emissions
            
        org_mass_pe += plus_mass_pe
        org_mass_se += plus_mass_se
            
        org_frac_se = calc_org_frac_se(org_mass_pe, org_mass_se)
            
        org_mass_se -= minus_mass *  org_frac_se
        org_mass_pe -= minus_mass *  (1-org_frac_se)
            
        dyn_res_pe[:,ts+1] = org_mass_pe
        dyn_res_se[:,ts+1] = org_mass_se
#    dyn_res_frac_se =  dyn_res_se /(dyn_res_se + dyn_res_pe)     
    return dyn_res_pe, dyn_res_se

def calc_org_frac_se(mass_pe, mass_se):
        X=[]
        for i in range(len(mass_pe)):
            if (mass_pe[i] + mass_se[i]) == 0:
                X.append(0)
            else:
                X.append(mass_se[i]/(mass_pe[i] + mass_se[i]))
        return asarray(X)
        
        

    

    
    
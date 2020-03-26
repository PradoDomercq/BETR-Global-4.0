# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:00:50 2013

@author: sebastian
"""
import pdb

def  x_to_cellcomp(x, cmap, ncomp ):
    x = cmap[x]
    cell = x/ncomp
    comp = x % ncomp
    return cell, comp    

def normalizefluxkey(flux_key):
    for i in range(len(flux_key)):
        for ii in range(len(flux_key[i])):
            for iii in range(len(flux_key[i][ii])):
                fsum=0.0
                for key in flux_key[i][ii][iii].keys():
                    fsum += flux_key[i][ii][iii][key]
                if fsum !=0.0:
                    for key in flux_key[i][ii][iii].keys():
                        flux_key[i][ii][iii][key] = flux_key[i][ii][iii][key]/fsum
    return flux_key
                    
    
    
def mk_airwater(m):
    raise NotImplementedError    
    
def mk_deg(m):
    raise NotImplementedError    
    
def mk_sed(m):
    raise NotImplementedError    
    
def mk_dep_air(m):
    raise NotImplementedError    
    
def mk_dep_water(m):
    raise NotImplementedError    
    
def mk_compmap(m):
    nr=0
    compmap = []
    for i in range(m.nocells*m.nocomp):
        if i in m.killidx:
            nr+=1
        else:
            compmap.append(nr)
            nr+=1
            
    return compmap 
    

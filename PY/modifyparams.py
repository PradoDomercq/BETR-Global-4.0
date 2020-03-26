################################################################################
'''Module "modifyparameters", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module calculates derived parameters, changes units etc.'''
################################################################################
from numpy import *
import os
import string
import inspect
import glob
import csv
from globalz import *
import sys

def _addarray(oldarray, plusarray):
    ''' attaches the contents of derived parameter array
    to the original parameter array''' 
    if oldarray.shape!= plusarray.shape:
        sys.exit('modparams: inconsistency in number of cells. Aborting!')
    par=empty(oldarray.shape, dtype=oldarray.dtype.descr+plusarray.dtype.descr)
    for nam in oldarray.dtype.names:
        par[nam]=oldarray[nam]
    for nam in plusarray.dtype.names:
        par[nam]=plusarray[nam]
    return(par)

def modparams(model):
    ''' Derive new parameters from original ones. Here, the scavenging
    ratio is selected with respect to snow / rain'''
    scavarr=(model.par['tair2'] < 273.15).astype(int)*model.par['scavsnow']\
             +(model.par['tair2'] >= 273.15).astype(int)*model.par['scavrain']
    scavarr.dtype=[('scavrat','f8')]
    par =_addarray(model.par, scavarr)
    return(par)
  

################################################################################
'''Module "readinput", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module reads and batters into shape the input parameters'''
################################################################################
from numpy import *
import os
import string
import inspect
import glob
import csv
from globalz import *
import sys
import pdb
import re

def _readEnvironment(partype, fn):
    """ Reads a file containing environmental variables.
    <partype> : 'const' or 'seasonal'. Determines whether the input file
    describes time-constant or time varying parameters. Both filetypes need
    to have a column 'CELL'. Seasonal-type parameter files need to have a
    column 'TS' (timestep).
    <filename> : the name of the parameter file.
    <envdir> : The directory of the parameter file.
    Returns record array (cells x 1) for partype='const' and
    (cells x timesteps) for partype='seasonal'.  
    """
    f=open(fn,'r')
    c='#'; line=''
    while c == "#":
        line_old=line
        line=f.readline()
        c=line[0]
    header=string.split(line_old)
    header[0]=header[0].lstrip("#")
    f.close()
    if partype == 'const':    
        dt=['int']+(len(header)-1)*['f8']
        pars=loadtxt(fn, dtype={'names':header, 'formats':dt}, comments='#')
        cells=unique(pars['CELL']).shape[0] #number of cells read
        ## remove field 'CELL' and sort according to increasing cell number
        pars=sort(pars,axis=0,order=['CELL'])
        fields=list(pars.dtype.names)
        fields.remove('CELL')
        pars=pars[fields]
        pars=reshape(pars,(cells,1))
        ## check dimension
        assert(cells == pars.shape[0])
    elif partype == 'seasonal':
        dt=2*['int']+(len(header)-2)*['f8']
        pars=loadtxt(fn, dtype={'names':header, 'formats':dt}, comments='#')
        cells=unique(pars['CELL']).shape[0]
        ts=unique(pars['TS']).shape[0]
        ## consistency check
        assert(cells*ts == pars.shape[0])
        ## remove field 'CELL' and sort according to
        ## 1) increasing cell number and 2) increasing timestep
        pars=sort(pars, axis=0, order=['CELL','TS'])
        fields=list(pars.dtype.names)
        fields=filter(lambda x: x!='CELL' and x!='TS', fields)
        pars=pars[fields]
        pars=reshape(pars,(cells, ts))
    return(pars)

def constructEnvironment(constfile, seasonalfile):
    '''Takes two filenames for constant and seasonally variable parameters
    respectively, and constructs one (cells x timesteps) record-array.
    This array contains all parameters that refer to single cells.
    Returns a tuple containing that array and a seasonally averaged array.
    '''
    parconst=_readEnvironment('const', constfile)
    parseasonal=_readEnvironment('seasonal', seasonalfile)

    ## check whether dimensions are compatible
    assert(parconst.shape[1] == 1)
    assert(parconst.shape[0] == parseasonal.shape[0])
    cells=parconst.shape[0]
    ts=parseasonal.shape[1]
    par=empty((cells,ts), dtype=parconst.dtype.descr+parseasonal.dtype.descr)
    for nam in parseasonal.dtype.names:
        par[nam]=parseasonal[nam]
    for nam in parconst.dtype.names:
        par[nam]=parconst[nam]
    return(par)

def readCompartments(fn):
    """ Reads the list of compartments considered by the model.
    Returns a dictionary containing their ID-numbers as keys and dictionaries
    {name: , temp_variable: [,other_column_names: ...]} as values."""
    compdict={}
    f=open(fn, 'r')
    header=[]
    while True:
        line=f.readline()
        if line[0] == "#":
            continue
        header=string.split(line)
        del(header[0])
        break
    while True:
        line=f.readline()
        if line == '': break
        if line[0] =='#' or line == '\n': continue
        line=string.split(line)
        num=line[0]
        del line[0]
        compdict[int(num)] = dict(zip(header,line))
    return(compdict)

def readFlows(complist, flowdir):
    #withoutQ
    """ Reads matrices of flows between cells.
    Reads all files in the directory <flowdir>, but uses only those,
    that define flows between compartments in <complist>.
    Returns a list that contains tuples (fromcompartment,tocompartment)
    followed by the respective flow-array in coordinate form."""
    fmatrixdict={}
    for fn in glob.iglob(os.path.join(flowdir,"*")):
        f=open(fn,'r')
        c='#'; line=''
        while c == "#":
            line_old=line
            line=f.readline()
            c=line[0]
        fro, to=line_old.split()
        fro=fro.lstrip("#")
        fro,to=[int(x) for x in[fro,to]]
        f.close()
        if not (fro in complist and to in complist):
            continue
        else:
            fmatrixdict[(fro,to)]=loadtxt(fn, comments='#')
    return(fmatrixdict)
    """
    SScehnker, we must somehow make sure that the minimum flush rate of freshwater 
    bodies is in the size of 5 if fresh water exists in a compartment. 
    
    Miller, J. R.;  Russell, G. L.;  Caliri, G. J. Climate 1994, 7, 914-928. 
    
    """

def readChemicals(chemlist, fn):
    ''' Reads the chemical properties of the chemical IDs in <chemlist>.
    Returns dictionary with IDs as keys and dictonaries with substance
    properties as values.
    <dbfile> is the file containing the database, and
    <chemdir> is its directory.'''
    f=open(fn,'r')
    # read header (fieldnames)
    c='#'; line=''; count=0
    while c == "#":
        line_old=line
        line=f.readline()
        c=line[0]
        count+=1
    header=string.split(line_old)
    header[0]=header[0].lstrip("#")
    f.seek(0)
    # jump over comments
    for i in range(0,count-1):
        f.readline()
    # construct dictionary
    r=csv.DictReader(f, fieldnames=header, delimiter=' ', quotechar='"')
    chemdict={}
    for chem in r:
        chemid=int(chem.pop('ID'))        
        if not chemid in chemlist:
            continue
        chemdict[chemid]=chem
        ## convert float-type fields, that is all except 'Name' and 'notes'
        for fname in header[2:25]:
            chemdict[chemid][fname]=float(chemdict[chemid][fname])
    f.close()
    return(chemdict)

def readProcesses(compdict, fn):
    proclist=[]
    f=open(fn,'r')
    count=0
    while True:
        line=f.readline()
        if line == '': break
        if line[0] == '#' or line == '\n': continue
        line=line.split()
        comps=[int(x) for x in line[1:]]
        if all([x in compdict.keys() for x in comps]):
            pl=[line[0]]
            pl.extend(comps)
            proclist.append(tuple(pl))
            count+=1
        ## check for double names of processes
        pnames=[x[0] for x in proclist]
        if len(unique(pnames)) != len(pnames):
            sys.exit('readinput.py: readProcesses: '
                  +'detected non-unique process names.\nAborting !')
    return(proclist)

def readControl(controlfile):
    controldict={}
    f=open(controlfile,'r')
    while True:
        line=f.readline()
        if line == '': break
        if line[0] =='#' or line == '\n': continue
        line=line.split()
        controldict[line[0]]=line[1]
    f.close()
    return(controldict)

    


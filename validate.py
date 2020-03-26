#!/usr/bin/python
################################################################################
'''Script "validate.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This script serves to validate BETR-Research against BETR-Global (VBA).
The validation is based on the system-matrices of D-values. The script has to
be executed from the directory where it resides. It takes the desired precision
and the IDs of the chemicals that should be used in the validation as arguments.
For example
validate.py 1e-13 2 11
will validate the model with chemicals 2 (PCB-8) and 11 (PCB-118) and fail if a
relative error, (Dvba-Dpy)/Dvba, greater than 1e-13 between any D-values from
BETR-Global (Dvba) and D-values from BETR-Research (Dpy) occurs.
The D-values from BETR-Global for each chemical have are contained in files
Book[n].txt, where [n] is the ID of the chemical. These files reside in the
subdirectoy "vbaresults".
'''
################################################################################
from numpy import *
import sys
from scipy import sparse
import cPickle
import subprocess
import os
import shutil
from pylab import *
import BETRS
reload(BETRS)
from BETRS import *
import helpers

# def cell2regcomp(cell):
#     reg=floor(cell/7) + 1
#     comp=mod(cell,7) + 1
#     return (reg.astype(int),comp.astype(int))

# def tocell(reg,comp):
#     cell=7*(reg-1)+comp-1 
#     return cell

def validate(level, chemicals):
    print('Validating BETR-Research for chemicals')
    print(chemicals)
    print('at level {0:.3e}'.format(level))
    valid=True
    difflist=[]
    for chem in chemicals:
        print('Calculating system matrices for chemical %d') % (chem)
        m=Model(chem,controlfile='control_validation.txt')
        shutil.copy(os.path.join('Output','default','matrixdump.cpk'),
                os.path.join('Validation'))
  
        ## load BETRS-output
        fn=os.path.join('Output','default','matrixdump.cpk')
        f=open(fn,'r')
        Dmatpy=cPickle.load(f)
        f.close()

        ## load BETR-VBA output
        fnvba=os.path.join('Validation','vbaresults','Book'+str(chem)+'.txt')
        try:
            Dmatvba=loadtxt(fnvba, dtype='S40',skiprows=2)
        except IOError:
            sys.exit('Did not find the file {0:s}. Aborting !'.format(fnvba))

        for mo in arange(0,12):
            print('comparing substance %s, month %i') % (chem,mo+1)
            Dres=Dmatpy[mo]
            ## get empty volume elements
            killidx=findemptycells(Dres)
            Dres=Dres.toarray()

            ## find irrelevant entries in VBA output
            keeplist = []
            for c in enumerate(Dmatvba[:,1]):
                try:
                    check = float(c[1])
                    keeplist.append(c[0])
                except ValueError:
                    pass
            Dmatvba=Dmatvba[array(keeplist),:].astype(float)
            ij=(Dmatvba[:,1].astype(int)-1, Dmatvba[:,0].astype(int)-1)
            Dvba=sparse.coo_matrix((Dmatvba[:,3+mo],ij),
                                   shape=(2016,2016)).toarray()
            ## fix erroneous non-zero entries on diagonal,
            ## check compatibility with killidx
            zerodiag = Dvba-diag(diag(Dvba))
            zerovols=where(sum(zerodiag,axis=0) == 0)[0]
            if all(zerovols == killidx):
                print('Test for identical zero-volumes: passed.')
            else:
                print('Failed test for identical zero volumes !')
                valid=False
            Dvba[zerovols,zerovols]=0
            Dvba=Dvba-diag(2*diag(Dvba))

            #### compare
            ## zero elements ?
            diff1=where(logical_and(Dres==0, Dvba!=0))
            diff2=where(logical_and(Dres!=0, Dvba==0))
            if diff1[0].size > 0 or diff2[0].size > 0:
                print('Validation failed ! There are zeros in one matrix that'\
                      +' are non-zero in the other:')
                print('research == 0, vba != 0: {0:s}'.format(diff1))
                print('research != 0, vba == 0: {0:s}'.format(diff2))
                valid=False
            else:
                print('Test for equality of sparsity pattern: passed.')
            ###########  OK ##########################
            di=(Dres-Dvba)/Dvba
            difflist.extend(list(di.ravel()[(logical_not(isnan(di.ravel())))]))
            di=abs(di)
            s1=where(di > level)
            fromc=list(cell2regcomp(s1[1],m)[1])
            toc=list(cell2regcomp(s1[0],m)[1])
            cells=list(cell2regcomp(s1[0],m)[0])
            ft=zip(fromc,toc)
            diffdict={}
            for i in enumerate(ft):
                if diffdict.has_key(i[1]):
                    diffdict[i[1]].append(cells[i[0]])
                else:
                    diffdict[i[1]]=[cells[i[0]]]
            if bool(diffdict):
                print('Validation not passed, relative differences > {0:.3e}'\
                      +' occured :'.format(level))
                print(diffdict)
                print('erroneous from-to:')
                print(unique(ft))
                valid=False
            else:
                print('Matices equal within rtol={0:.3e}'.format(level))
    print(60*'*')
    if valid:
        print('Overall validation: PASSED')
    else:
        print('Overall validation: FAILED')
        
    flatdiffs=array(difflist)
    return(flatdiffs)

if __name__ == '__main__':
    level=float(sys.argv[1])
    chemicals=[int(x) for x in sys.argv[2:]]
    diffs=validate(level, chemicals)
    hist(diffs,30)
    show()

################################################################################

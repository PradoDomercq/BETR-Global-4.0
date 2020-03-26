################################################################################
'''Module "output.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module contains output-routines for BETR-Research'''
################################################################################
from numpy import *
import sys
import os
import csv
import cPickle
import pdb
#import write_ncfile
#reload(write_ncfile)
#try:
from write_ncfile import *
#except ImportError:
#	from write_ncfile_old import *

from helpers import *
from helpers3 import *
from secondary_emission_reconstruction import reconstruct_dyn_res

def write_cpk_file(fn, out):
    
    fncpk=fn+'_out.cpk'
    if not os.path.exists(os.path.dirname(fncpk)):
        os.mkdir(os.path.dirname(fncpk))
    else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
    f=open(fncpk,'w')
    cPickle.dump(out,f)
    f.close()

def write_output_summary(m, fn, nchemical, nrun, nchemdb, nseasonalparfile, nconstantparfile,\
                       ncompartmentfile, nflowdir, nprocessfile, ncontrolfile):
    if not os.path.exists(os.path.dirname(fn)):
        os.mkdir(os.path.dirname(fn))    
    writer = csv.writer(open(fn, 'w'), delimiter = ' ')
    writer.writerows([["BETR-Global 3.0"],
                      ['runID', nrun], 
                      ['seasfile', nseasonalparfile],
                      ['chemdata', nchemdb], 
                      ['chemnr', nchemical],
                      ['compfile', ncompartmentfile], 
                      ['constparfile', nconstantparfile],
                      ['flowdirectory', nflowdir], 
                      ['procfile', nprocessfile],
                      ['contfile', ncontrolfile]])      
      
def write_output_summary2(m, fn, nemfile):
    writer = csv.writer(open(fn, 'a'), delimiter = ' ')
    writer.writerow(['emisfile', nemfile])      
    
def write_output_summary3(m, fn, nsolvfile):
    writer = csv.writer(open(fn, 'a'), delimiter = ' ')
    writer.writerow(['solvfile', nsolvfile])    
    
def write_output_ss(m, fn, units, netcdf, cpk):     # parameter cpk added by HW
    ''' Output for steady-state results'''
    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero-volumes / Z-values
        varr=mean(m.vdict[c]['bulk'],axis=1)
        zv_arr=mean(m.zdict[c]['bulk'],axis=1)*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        if 'mol' in units:
            out[c]['mol']=m.ss_res[idx]
        if 'kg' in units:
            out[c]['kg']=m.ss_res[idx]*m.chemdict['molmass']/1000.0
        if 'mol_per_m3' in units:
            out[c]['mol_per_m3']=m.ss_res[idx]/varr
            out[c]['mol_per_m3'][zerovolidx]=0
        if 'kg_per_m3' in units:
            out[c]['kg_per_m3']=m.ss_res[idx]*m.chemdict['molmass']/1000.0/varr
            #out[c]['pg_per_m3']=m.ss_res[idx]/varr/8760*1e15   #25.08.2016 fy
            out[c]['kg_per_m3'][zerovolidx]=0
        if 'Pa' in units:
            out[c]['Pa']=m.ss_res[idx]/zv_arr
            out[c]['Pa'][zerozvidx]=0
            
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):         # moved upwards
        os.mkdir(os.path.dirname(fn))
            
    if cpk:                                             # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            #outarray=zeros((m.nocomp,nlon,nlat)) #12,24
            outarray=zeros((m.nocomp,nlat,nlon)) #12,24 # Modified by RKG, 28.07.2014
            for c in m.compdict.keys():
                #outarray[c-1,:,:]=out[c][u].reshape(nlon,nlat) #12,24
                outarray[c-1,:,:]=out[c][u].reshape(nlat,nlon) #12,24 # Modified by RKG, 28.07.2014
            writenc(outarray,fnnc,False,True,1,varname='V',unit=u)

def write_output_dyn(m, fn, units, netcdf, cpk):     # parameter cpk added by HW
    ''' Output for dynamic results'''
    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    #print 'write_output_dyn: self.dyn_res:', type(m.dyn_res),m.dyn_res.shape # RKG, 03.07.2014
    timesteps=m.dyn_res.shape[1]
    #print 'write_output_dyn: timesteps = ', timesteps # RKG, 03.07.2014
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)), 
                       tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99

        if 'mol' in units:
            out[c]['mol']=m.dyn_res[idx]
        if 'kg' in units:
            out[c]['kg']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0
        if 'mol_per_m3' in units:
            out[c]['mol_per_m3']=m.dyn_res[idx]/varr
            out[c]['mol_per_m3'][zerovolidx]=0
        #if 'kg_per_m3' in units:
        if  'ng_per_m3' in units:
            #out[c]['kg_per_m3']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0/varr
            #out[c]['kg_per_m3'][zerovolidx]=0
            out[c]['ng_per_m3']=m.dyn_res[idx]*m.chemdict['molmass']*1e9/varr  #for EST_45_3349_2011   FIG.4(b)
            out[c]['ng_per_m3'][zerovolidx]=0
        if 'Pa' in units:
            out[c]['Pa']=m.dyn_res[idx]/zv_arr
            out[c]['Pa'][zerozvidx]=0
    
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):   # moved upwards
        os.mkdir(os.path.dirname(fn))
        
    if cpk:                          # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,m.nocomp,nlat,nlon))# 12,24
            for c in m.compdict.keys():
                for t in arange(0,timesteps):
                    outarray[t,c-1,:,:]=out[c][u][:,t].reshape(nlat,nlon) #12,24
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)


def write_output_se(m, fn, units, netcdf, cpk, scenario="air"):     # parameter cpk added by HW
    ''' Output for secondary emissions results'''
#    from secondary_emission_reconstruction import reconstruct_dyn_res
    units = ["pe", "se", 'pe_frac', 'se_frac']
    dyn_res_pe, dyn_res_se =  reconstruct_dyn_res(m, scenario=scenario,nseasons=12)
    
    frac_se = dyn_res_se / (dyn_res_pe + dyn_res_se)
    frac_pe = dyn_res_pe / (dyn_res_pe + dyn_res_se)
    
    index_NaN = isnan(frac_se)
    frac_se[index_NaN] = 0
    frac_pe[index_NaN] = 0

    nlat=int(sqrt(m.matdim/m.nocomp/2))                  # added by HW
    nlon=int(sqrt(m.matdim/m.nocomp/2)*2)                # added by HW
    out={}
    timesteps=m.dyn_res.shape[1]
    #print 'write_output_se: timesteps = ', timesteps # RKG, 03.07.2014
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    for c in m.compdict.keys():
        out[c]={}
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)), 
                       tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out[c]['pe']=dyn_res_pe[idx]
        out[c]['se']=dyn_res_se[idx]
        out[c]['pe_frac']=frac_pe[idx]
        out[c]['se_frac']=frac_se[idx]
        

#        if 'mol' in units:
#            out[c]['mol']=m.dyn_res[idx]
#        if 'kg' in units:
#            out[c]['kg']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0
#        if 'mol_per_m3' in units:
#            out[c]['mol_per_m3']=m.dyn_res[idx]/varr
#            out[c]['mol_per_m3'][zerovolidx]=0
#        if 'kg_per_m3' in units:
#            out[c]['kg_per_m3']=m.dyn_res[idx]*m.chemdict['molmass']/1000.0/varr
#            out[c]['kg_per_m3'][zerovolidx]=0
#        if 'Pa' in units:
#            out[c]['Pa']=m.dyn_res[idx]/zv_arr
#            out[c]['Pa'][zerozvidx]=0
    
    if not out:
        print("WARNING: empty output requested !\n")
        
    if not os.path.exists(os.path.dirname(fn)):   # moved upwards
        os.mkdir(os.path.dirname(fn))
        
    if cpk:                          # condition added by HW
        fncpk=fn+'_out.cpk'
        #if not os.path.exists(os.path.dirname(fncpk)):
            #os.mkdir(os.path.dirname(fncpk))
        #else:
        if os.path.exists(fncpk):
            print('Attention: overwriting %s\n') % (fncpk)
        f=open(fncpk,'w')
        cPickle.dump(out,f)
        f.close()
    if netcdf:
        for u in units:
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,m.nocomp,nlat,nlon))# 12,24
            for c in m.compdict.keys():
                for t in arange(0,timesteps):
                    outarray[t,c-1,:,:]=out[c][u][:,t].reshape(nlat,nlon) #12,24
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)
            
            
def write_output_dyn_tmp(m, fn, timesteps):                             # function added by HW
    nreg=m.matdim/m.nocomp    
    out_tmp = zeros([nreg*m.nocomp, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), m.nocomp)#289
    out_tmp[:, 1] = repeat(range(1, m.nocomp+1), nreg)#288
    
    #timesteps=m.dyn_res.shape[1]
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)),
                        tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.dyn_res[idx, timesteps-1] # unit = mol   #288 
          
    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')  # changed from .8e to .16e
    
    
def write_output_end_txt(m, fn, timesteps):                             # function added by HW
    nreg=m.matdim/m.nocomp
    out_tmp = zeros([nreg*m.nocomp, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), m.nocomp)#289
    out_tmp[:, 1] = repeat(range(1, m.nocomp+1), nreg)#288
    
    #timesteps = m.dyn_res.shape[1]    
    periods=int((timesteps-1)/float(m.nots))
    if periods != (timesteps-1)/float(m.nots):
        sys.exit('timesteps of output {0:d} not multiple of '
                 +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
        varr=hstack((ones((idx.shape[0],1)),
                     tile(m.vdict[c]['bulk'], (1,periods))))
        zv_arr=hstack((ones((idx.shape[0],1)),
                        tile(m.zdict[c]['bulk'], (1,periods))))*varr
        zerovolidx=where(varr==0)
        varr[zerovolidx]=-9999.99
        zerozvidx=where(zv_arr==0)
        zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.dyn_res[idx, timesteps-1] # unit = mol   #288  #
          
    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')   # changed from .8e to .16e

def write_output_ss_txt(m, fn):                             # function added by RKG, 04.01.2014
    nreg=m.matdim/m.nocomp
    out_tmp = zeros([nreg*m.nocomp, 3])#288*7
    out_tmp[:, 0] = tile(range(1, nreg+1), m.nocomp)#289
    out_tmp[:, 1] = repeat(range(1, m.nocomp+1), nreg)#288
    
    #timesteps = m.dyn_res.shape[1]    
    #periods=int((timesteps-1)/float(m.nots))
    #if periods != (timesteps-1)/float(m.nots):
    #    sys.exit('timesteps of output {0:d} not multiple of '
    #             +'seasons ({1:d}): Aborting !'.format(timesteps,m.nots))
    
    for c in m.compdict.keys():
        idx=arange(c-1,m.matdim,len(m.compdict.keys()))
        ## deal with zero volumes / Z-values
    #    varr=hstack((ones((idx.shape[0],1)),
    #                 tile(m.vdict[c]['bulk'], (1,periods))))
    #    zv_arr=hstack((ones((idx.shape[0],1)),
    #                    tile(m.zdict[c]['bulk'], (1,periods))))*varr
    #    zerovolidx=where(varr==0)
    #    varr[zerovolidx]=-9999.99
    #    zerozvidx=where(zv_arr==0)
    #    zv_arr[zerozvidx]=-9999.99
        
        out_tmp[(c-1)*nreg : (c-1)*nreg+nreg, 2] = m.ss_res[idx] # unit = mol   #288  #

    savetxt(fn, out_tmp, fmt = ['%i', '%i', '%.16e'], delimiter=' ')   

    
def write_bigDlist_txt(bigDlist, fn2, fn3, m):             # function added by HW
    nreg=m.matdim/m.nocomp
    #dtypeDlist = [('regfrom', int), ('compfrom', int), ('regto', int), ('compto', int), 
                      #(str(1), float), (str(2), float), (str(3), float), (str(4), float), (str(5), float), (str(6), float),
                      #(str(7), float), (str(8), float), (str(9), float), (str(10), float), (str(11), float), (str(12), float)]
    Dflow = zeros((7*nreg*nreg, 16))
    fromreg = hstack([sort(range(1,nreg+1)*nreg)]*7)
    toreg = hstack([range(1,nreg+1)*nreg]*7)
    Dflow[:,range(0,4)] = transpose([fromreg,
                                     hstack([[1]*nreg*nreg, [2]*nreg*nreg, [4]*nreg*nreg, [5]*nreg*nreg, [7]*nreg*nreg, [8]*nreg*nreg, [9]*nreg*nreg]),
                                     toreg,
                                     hstack([[1]*nreg*nreg, [2]*nreg*nreg, [4]*nreg*nreg, [5]*nreg*nreg, [7]*nreg*nreg, [8]*nreg*nreg, [9]*nreg*nreg])])
    Dflow = Dflow[Dflow[:,0] != Dflow[:,2], :]
    for month in range(0,12):
        Dflow[:,month + 4] = bigDlist[month][tocell(Dflow[:,2], Dflow[:,3], m), tocell(Dflow[:,0], Dflow[:,1], m)]
    Dflow = Dflow[findnonemptycells(Dflow[:, range(4,16)]), :]
    
    Ddiag = zeros((m.nocomp*nreg, 16))
    fromreg = sort(range(1,nreg+1)*m.nocomp)
    fromcomp = range(1,m.nocomp+1)*nreg
    Ddiag[:,range(0,4)] = transpose([fromreg, fromcomp, fromreg, fromcomp])
    for month in range(0,12):
        Ddiag[:,month + 4] = bigDlist[month][tocell(Ddiag[:,2], Ddiag[:,3], m), tocell(Ddiag[:,0], Ddiag[:,1], m)]
    Ddiag = Ddiag[findnonemptycells(Ddiag[:, range(4,16)]), :]
    
    Dint = zeros((nreg*m.nocomp*m.nocomp, 16))
    fromreg = sort(range(1,nreg+1)*m.nocomp*m.nocomp)
    fromcomp = hstack([sort(range(1,m.nocomp+1)*m.nocomp)]*nreg)
    tocomp = range(1,m.nocomp+1)*m.nocomp*nreg
    Dint[:,range(0,4)] = transpose([fromreg, fromcomp, fromreg, tocomp])
    Dint = Dint[Dint[:,1] != Dint[:,3],:]
    Dint = Dint[Dint[:,0] == Dint[:,2],:]   # unneccesary
    for month in range(0,12):
        Dint[:,month + 4] = bigDlist[month][tocell(Dint[:,2], Dint[:,3], m), tocell(Dint[:,0], Dint[:,1], m)]
    Dint = Dint[findnonemptycells(Dint[:, range(4,16)]), :]
    
    bigDlist2 = vstack([Dflow, Ddiag, Dint])
    '''
    for i in range(len(bigDlist2)):
        if bigDlist2[i][0]==1 and bigDlist2[i][1]==1:
            print bigDlist2[i][0:5]
        if bigDlist2[i][2]==1 and bigDlist2[i][3]==1:
            print bigDlist2[i][0:5]
    pdb.set_trace()
    '''
    
    #output timelist FY
    timelist = zeros((len(bigDlist2),16))    
    for i in range(len(bigDlist2)):
        zfromcell=bigDlist2[i][0]-1;zfromcomp=bigDlist2[i][1]

        timelist[i][0:4]=bigDlist2[i][0:4]
        vz=m.vdict[zfromcomp]['bulk'][zfromcell]*m.zdict[zfromcomp]['bulk'][zfromcell]
        for j in range(12):
            if bigDlist2[i][4+j]==0:
                timelist[i][4+j]=0
            else:
                timelist[i][4+j]=abs(vz[j]/bigDlist2[i][4+j])

    '''
    for i in range(len(timelist)):
        if timelist[i][0]==1 and timelist[i][1]==1:
            print timelist[i][0:5]
        if timelist[i][2]==1 and timelist[i][3]==1:
            print timelist[i][0:5]
    pdb.set_trace()
    '''

    savetxt(fn2, bigDlist2, 
            fmt = ['%i', '%i', '%i', '%i', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e'],            
            delimiter=' ')
    savetxt(fn3, timelist, 
            fmt = ['%i', '%i', '%i', '%i', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e', '%.8e'],            
            delimiter=' ')


def write_output_dflux(m,fn,netcdf=False,units=['mol'],nlat=12,nlon=24):
    '''Output all files related to the integrated fluxes'''                                
#    possiblefiles = ["flow", "airwater", "deg", "sed", "dep_air", "dep_water"]
    
    out = mk_allfluxes(m)
    k = out.keys()[0]
    timesteps=len(out[k])
    write_cpk_file(fn, out)
    if netcdf:
        for u in out.keys():
            fnnc=fn+'_'+u+'.nc'
            if os.path.exists(fnnc):
                print('Attention: overwriting %s\n') % (fnnc)
            outarray=zeros((timesteps,1,nlat,nlon))
            for t in arange(0,timesteps):
                    outarray[t,0,:,:]=out[u][t].reshape(nlat,nlon)
            writenc(outarray, fnnc, True, True, 1, varname='V',unit=u)
##############################################################################3
#SSchenker    The following functions are needed for Linear Approx. Fluxes
#             and for the follow_flow mode (Correct fluxes)
    
def mk_allfluxes(m, nseasons=12):
    '''Calculates the net in/outflow from all compartments at every timestep'''
    emiss=get_emission_array(m) #read emission which will be used in checking the mass balance
    flux_key = normalizefluxkey(m.flux_key)
    # Calculate the inventory in each compartment
    inventory = zeros((13,m.nocells,m.nocomp))
    for c in range(m.nocomp):
        idx=arange(c,m.matdim,m.nocomp)
        for cell in range(len(idx)):
            for ts in range(13):            
                inventory[ts][cell][c]=m.dyn_res[idx[cell]][ts]
    
    fluxdict = {}
    fluxes=[
            "flow_in_upper_gas", "flow_out_upper_gas","upper_gas_sed","upper_gas_deg", 
            "upper_gas_to_all", "all_to_upper_gas","upper_gas", "upper_gas_ratio",
            # "upper_gas_ratio"=input/output, which is used to check the mass balance equation. If the ratio is 1, it means the input equals the output. 
            "flow_in_lower_gas", "flow_out_lower_gas","lower_gas_sed","lower_gas_deg", 
            "lower_gas_to_all", "all_to_lower_gas","lower_gas", "lower_gas_ratio",
            "flow_in_upper_faerosol","flow_out_upper_faerosol","flow_in_upper_caerosol","flow_out_upper_caerosol",
            "upper_caerosol_sed","upper_caerosol_deg",
            "upper_caerosol_to_all","all_to_upper_caerosol",
            "upper_caerosol","upper_caerosol_ratio",
            "flow_in_lower_faerosol","flow_out_lower_faerosol","flow_in_lower_caerosol","flow_out_lower_caerosol",
            "lower_faerosol_sed","lower_caerosol_sed","lower_faerosol_deg","lower_caerosol_deg",
            "lower_faerosol_to_all","all_to_lower_faerosol","lower_caerosol_to_all","all_to_lower_caerosol",
            "lower_faerosol","lower_caerosol","lower_faerosol_ratio","lower_caerosol_ratio",
            ]

    for key in fluxes:
        fluxdict[key] = zeros((len(m.flux_res), m.nocells))   #len(m.flux_res)=12    
    compmap = mk_compmap(m)
    ts0 = time.time()
    ts1 = time.time()
    for ts in range(len(m.flux_res)):
        season = ts % nseasons
        year = ts / nseasons
        if season ==0:
            print ("Processing year %i \t" %(year)),
        
        fluxmat = m.flux_res[ts]
        x,y = fluxmat.nonzero()
        for i in range(len(x)):                
            value = fluxmat[x[i],y[i]]            
            xcell, xcomp = x_to_cellcomp( x[i], compmap, m.nocomp)
            ycell, ycomp = x_to_cellcomp( y[i], compmap, m.nocomp)
            #print 'xcell, xcomp',xcell, xcomp
            #print 'value', value
            
            if x[i]==y[i]:
                try :
                    match_dict = m.flux_key[season][xcell][xcomp]
                    #print 'match_dict', match_dict
                except KeyError:
                    print [[season],[xcell],[xcomp]]
                for mnkey in match_dict.keys():
                    
                    if type(mnkey) is tuple:                             
                        if xcomp in [0]:
                            fluxdict["flow_in_upper_gas"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_upper_gas"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value                            
                        if xcomp in [1]:
                            fluxdict["flow_in_lower_gas"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_lower_gas"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [7]:
                            fluxdict["flow_in_upper_caerosol"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_upper_caerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [8]:
                            fluxdict["flow_in_lower_faerosol"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_lower_faerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp in [9]:
                            fluxdict["flow_in_lower_caerosol"][ts][mnkey[1]-1] -= flux_key[season][xcell][xcomp][mnkey] * value
                            fluxdict["flow_out_lower_caerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value 
                
                    if mnkey  == "sed":
                        if xcomp == 0:
                            fluxdict["upper_gas_sed"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 1:
                            fluxdict["lower_gas_sed"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 7:
                            fluxdict["upper_caerosol_sed"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 8:
                            fluxdict["lower_faerosol_sed"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 9:
                            fluxdict["lower_caerosol_sed"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                            
                    if mnkey  == "deg":
                        if xcomp == 0:
                            fluxdict["upper_gas_deg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 1:
                            fluxdict["lower_gas_deg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 7:
                            fluxdict["upper_caerosol_deg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 8:
                            fluxdict["lower_faerosol_deg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                        if xcomp == 9:
                            fluxdict["lower_caerosol_deg"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                            
                    if xcomp in [0] and mnkey in [1,2,3,4,5,6,7,8,9]:
                        fluxdict["upper_gas_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [1,2,3,4,5,6,7,8,9] and mnkey in [0]:
                        fluxdict["all_to_upper_gas"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [1] and mnkey in [0,2,3,4,5,6,7,8,9]:
                        fluxdict["lower_gas_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [0,2,3,4,5,6,7,8,9] and mnkey in [1]:
                        fluxdict["all_to_lower_gas"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [7] and mnkey in [0,1,2,3,4,5,6,8,9]:
                        fluxdict["upper_caerosol_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [0,1,2,3,4,5,6,8,9] and mnkey in [7]:
                        fluxdict["all_to_upper_caerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [8] and mnkey in [0,1,2,3,4,5,6,7,9]:
                        fluxdict["lower_faerosol_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [0,1,2,3,4,5,6,7,9] and mnkey in [8]:
                        fluxdict["all_to_lower_faerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [9] and mnkey in [0,1,2,3,4,5,6,7,8]:
                        fluxdict["lower_caerosol_to_all"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value
                    if xcomp in [0,1,2,3,4,5,6,7,8] and mnkey in [9]:
                        fluxdict["all_to_lower_caerosol"][ts][xcell] -= flux_key[season][xcell][xcomp][mnkey] * value

                    

        for i in range(m.nocells):
            
            fluxdict["upper_gas"][ts][i]=(fluxdict["flow_in_upper_gas"][ts][i]+fluxdict["all_to_upper_gas"][ts][i]+(emiss[ts][i][0]*730))\
                                          -(inventory[(ts+1)][i][0]-inventory[(ts)][i][0])\
                                          -(fluxdict["flow_out_upper_gas"][ts][i]+fluxdict["upper_gas_to_all"][ts][i]\
                                            +fluxdict["upper_gas_sed"][ts][i]+fluxdict["upper_gas_deg"][ts][i])
            fluxdict["upper_gas_ratio"][ts][i]=(fluxdict["flow_in_upper_gas"][ts][i]+fluxdict["all_to_upper_gas"][ts][i]+(emiss[ts][i][0]*730)\
                                          -(inventory[(ts+1)][i][0]-inventory[(ts)][i][0]))\
                                          /(fluxdict["flow_out_upper_gas"][ts][i]+fluxdict["upper_gas_to_all"][ts][i]\
                                            +fluxdict["upper_gas_sed"][ts][i]+fluxdict["upper_gas_deg"][ts][i])
            fluxdict["lower_gas"][ts][i]=(fluxdict["flow_in_lower_gas"][ts][i]+fluxdict["all_to_lower_gas"][ts][i]+(emiss[ts][i][1]*730))\
                                          -(inventory[(ts+1)][i][1]-inventory[ts][i][1])\
                                          -(fluxdict["flow_out_lower_gas"][ts][i]+fluxdict["lower_gas_to_all"][ts][i]\
                                            +fluxdict["lower_gas_sed"][ts][i]+fluxdict["lower_gas_deg"][ts][i])
            fluxdict["lower_gas_ratio"][ts][i]=(fluxdict["flow_in_lower_gas"][ts][i]+fluxdict["all_to_lower_gas"][ts][i]+(emiss[ts][i][1]*730)\
                                          -(inventory[(ts+1)][i][1]-inventory[ts][i][1]))\
                                          /(fluxdict["flow_out_lower_gas"][ts][i]+fluxdict["lower_gas_to_all"][ts][i]\
                                            +fluxdict["lower_gas_sed"][ts][i]+fluxdict["lower_gas_deg"][ts][i])
            fluxdict["upper_caerosol"][ts][i]=(fluxdict["flow_in_upper_caerosol"][ts][i]+fluxdict["all_to_upper_caerosol"][ts][i]+(emiss[ts][i][7]*730))\
                                          -(inventory[(ts+1)][i][7]-inventory[(ts)][i][7])\
                                          -(fluxdict["flow_out_upper_caerosol"][ts][i]+fluxdict["upper_caerosol_to_all"][ts][i]\
                                            +fluxdict["upper_caerosol_sed"][ts][i]+fluxdict["upper_caerosol_deg"][ts][i])
            fluxdict["upper_caerosol_ratio"][ts][i]=(fluxdict["flow_in_upper_caerosol"][ts][i]+fluxdict["all_to_upper_caerosol"][ts][i]+(emiss[ts][i][7]*730)\
                                          -(inventory[(ts+1)][i][7]-inventory[(ts)][i][7]))\
                                          /(fluxdict["flow_out_upper_caerosol"][ts][i]+fluxdict["upper_caerosol_to_all"][ts][i]\
                                            +fluxdict["upper_caerosol_sed"][ts][i]+fluxdict["upper_caerosol_deg"][ts][i])
            fluxdict["lower_faerosol"][ts][i]=(fluxdict["flow_in_lower_faerosol"][ts][i]+fluxdict["all_to_lower_faerosol"][ts][i]+(emiss[ts][i][8]*730))\
                                          -(inventory[(ts+1)][i][8]-inventory[(ts)][i][8])\
                                          -(fluxdict["flow_out_lower_faerosol"][ts][i]+fluxdict["lower_faerosol_to_all"][ts][i]\
                                            +fluxdict["lower_faerosol_sed"][ts][i]+fluxdict["lower_faerosol_deg"][ts][i])
            fluxdict["lower_faerosol_ratio"][ts][i]=(fluxdict["flow_in_lower_faerosol"][ts][i]+fluxdict["all_to_lower_faerosol"][ts][i]+(emiss[ts][i][8]*730)\
                                          -(inventory[(ts+1)][i][8]-inventory[(ts)][i][8]))\
                                          /(fluxdict["flow_out_lower_faerosol"][ts][i]+fluxdict["lower_faerosol_to_all"][ts][i]\
                                            +fluxdict["lower_faerosol_sed"][ts][i]+fluxdict["lower_faerosol_deg"][ts][i])
            fluxdict["lower_caerosol"][ts][i]=(fluxdict["flow_in_lower_caerosol"][ts][i]+fluxdict["all_to_lower_caerosol"][ts][i]+(emiss[ts][i][9]*730))\
                                          -(inventory[(ts+1)][i][9]-inventory[(ts)][i][9])\
                                          -(fluxdict["flow_out_lower_caerosol"][ts][i]+fluxdict["lower_caerosol_to_all"][ts][i]\
                                            +fluxdict["lower_caerosol_sed"][ts][i]+fluxdict["lower_caerosol_deg"][ts][i])
            fluxdict["lower_caerosol_ratio"][ts][i]=(fluxdict["flow_in_lower_caerosol"][ts][i]+fluxdict["all_to_lower_caerosol"][ts][i]+(emiss[ts][i][9]*730)\
                                          -(inventory[(ts+1)][i][9]-inventory[(ts)][i][9]))\
                                          /(fluxdict["flow_out_lower_caerosol"][ts][i]+fluxdict["lower_caerosol_to_all"][ts][i]\
                                            +fluxdict["lower_caerosol_sed"][ts][i]+fluxdict["lower_caerosol_deg"][ts][i])
                    
        print ("."),
        if season == (nseasons-1):
            print('  [%.3f s]')  % (time.time()-ts1)
            ts1=time.time()        
    print ("Time in mk_allfluxes(): %.3f min " % ((time.time()-ts0)/60.))
    return fluxdict        
                      
                          
#################################################################################
def get_emission_array(m):                #FY
    fn=os.path.join('Emissions', m.emissionfile)
    f=open(fn, 'r')
    lines=f.readlines()
    f.close

    emiss=zeros((12,m.nocells,m.nocomp))
    for lin in lines:
            if (lin[0] == '#' or lin == ''):
                continue
            [m,r,c1,c2,c3,val,perc1,perc2]=[x.lstrip() for x in \
                         [x.rstrip() for x in lin.split()]]
            
            emiss[int(m)-1][int(r)-1][int(c1)-1]=float(val)*float(perc1)
            emiss[int(m)-1][int(r)-1][int(c2)-1]=float(val)*(1-float(perc1))*float(perc2)
            emiss[int(m)-1][int(r)-1][int(c3)-1]=float(val)*(1-float(perc1))*(1-float(perc2))
    return emiss
        

        

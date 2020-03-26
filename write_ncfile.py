################################################################################
'''Module "write_ncfile.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>


module to write an .nc - file
input: lat : array of latitudes, lon : array of longitudes
        hastime: there is a time-dimension
        hascompartments: there is a compartment dimension
        monthspertimestep: months per timestep
        d = array of shape ([size(t),] [size(comps),] size(lat), size(lon))
        filename of nc-file
        varname: nc-name of variable (not "data")
        unit: unit of variable'''
################################################################################

from numpy import *
from netCDF4 import Dataset
import time
import sys
import pdb

def writenc(d, filename, hastime, hascompartments, monthspertimestep,
            #lat=arange(82.5,-90.0,-15)
            #lon=arange(-172.5, 180.0,15)
            varname='v', unit='nounit'):
    lat=arange(90.0-(180.0/d.shape[len(d.shape)-2]/2.0), -90.0, -180.0/d.shape[len(d.shape)-2])
    lon=arange(-180.0+(360.0/d.shape[len(d.shape)-1]/2.0), 180.0, 360.0/d.shape[len(d.shape)-1])
    #f=Dataset(filename, 'w', format='NETCDF3_64BIT')
    #f=Dataset(filename, 'w', format='NETCDF4')
    f=Dataset(filename,
              'w',
              #format='NETCDF3_CLASSIC',
              )
    f.Conventions='CF-1.1'
    f.History= 'Created '+time.ctime(time.time())
    
    latbnds = zeros((len(lat),2))
    lonbnds = zeros((len(lon),2))
    latbnds[:,0] = lat - 0.5*180/len(lat)
    latbnds[:,1] = lat +  0.5*180/len(lat) 
    lonbnds[:,0] = lon - 0.5*360/len(lon)
    lonbnds[:,1] = lon +  0.5*360/len(lon)  

    #sanity check lat/lon dimensions
    if (size(lat) != d.shape[len(d.shape)-2]) \
           or (size(lon) != d.shape[len(d.shape)-1]):
        sys.exit("inconsistent dimensions. aborting !")
    
    nextdim=0
    dimensions=[]
    if hastime:
        timedim=f.createDimension('time',d.shape[nextdim])
        timevar=f.createVariable('time','i4',dimensions=('time',))
        timevar.units='months'
        timevar.axis='T'
        timevar[:]=arange(0,d.shape[nextdim]*monthspertimestep, dtype='i4')
        dimensions.append('time')
        nextdim+=1

    if hascompartments:
        compdim=f.createDimension('compartment',d.shape[nextdim])
        compvar=f.createVariable('compartment','i4',
                                 dimensions=('compartment',))
        compvar.axis='Z'
        compvar[:]=arange(1,d.shape[nextdim]+1,dtype='i4')
        dimensions.append('compartment')
        nextdim+=1
 
        
    latdim=f.createDimension('lat',d.shape[nextdim])
    latvar=f.createVariable('lat','f8',('lat',))
    latvar.units='degrees_north'
    latvar.long_name='latitude'
    latvar.axis='Y'
    latvar.bounds='lat_bnds'
    latvar[:]=lat
    dimensions.append('lat')
    nextdim+=1

    londim=f.createDimension('lon',d.shape[nextdim])
    lonvar=f.createVariable('lon','f8',('lon',))
    lonvar.units='degrees_east'
    lonvar.long_name='longitude'
    lonvar.axis='X'
    lonvar.bounds='lon_bnds'
    lonvar[:]=lon
    dimensions.append('lon')

    nvdim=f.createDimension('nv',2)
    nvvar=f.createVariable('nv','i4',('nv',))
    latbndvar=f.createVariable('lat_bnds','f8',('lat','nv'))
    latbndvar[:]=latbnds
    lonbndvar=f.createVariable('lon_bnds','f8',('lon','nv'))
    lonbndvar[:]=lonbnds
    
 
    variable=f.createVariable(varname,'f8',dimensions)
    variable.units=unit
    variable[:]=d
    
    f.sync()
    f.close()

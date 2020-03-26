################################################################################
#    BETR-Research
#
#    A framework to create spatially resolved multimedia fate- and transport
#    models for chemical contaminants.
#
#    Copyright (C) 2010  Harald von Waldow <hvwaldow@chem.ethz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

from numpy import *
import sys
import pdb
import cPickle
import os



## these imports reload modules in each run, for making development
## workflow smoother
## adding the subdirectory "PY" to the module search path. "PY" contains
## model-specific code.

pydir=os.path.abspath(os.path.join(os.path.dirname(__file__),'PY'))
if not pydir in sys.path:
    sys.path.insert(0,pydir)

import readinput
import globalz
import zvalues
import tempcorrect
import volumes
import processes
import mkflowD
import mkflowDQ
import mkDmixQ
import mkbigD
import emissions
import helpers
import solver
import output
import modifyparams
reload(readinput)
reload(globalz)
reload(zvalues)
reload(tempcorrect)
reload(volumes)
reload(processes)
reload(mkflowD)
reload(mkflowDQ)
reload(mkDmixQ)
reload(mkbigD)
reload(emissions)
reload(helpers)
reload(solver)
reload(output)
reload(modifyparams)
from readinput import *
from globalz import *
from zvalues import *
from tempcorrect import *
from volumes import *
from processes import *
from mkflowD import *
from mkflowDQ import *
from mkDmixQ import *
from mkbigD import *
from emissions import *
from helpers import *
from solver import *
from output import *
from modifyparams import *
from LinFluxApprox import * 
import scipy.sparse as sp
##########################

class Model():
    """
    The class :class:`Model` is the user\'s interface to BETR-Research.
    
    :param int chemical: The ID of the chemical to be modelled (according to :file:`Chemicals/{chemdb}`)
    :param string run: String describing a particular model run. Used for the construction of output file paths.
    :param string chemdb: Filename of the chemicals database.
    :param string seasonalparfile: Filename of the parameter-file for seasonally changing values.
    :param string constantparfile: Filename of parameter-file for time-constant parameters.
    :param string compartmentfile: Filename of file describing the compartment.
    :param string flowdir: Name of the directory containing the convective flow descriptions.
    :param string processfile: Filename of the list of processes to be considered.
    :param string controlfile: Filename of the file containing options for this model run.

    """
    
    def __init__(self,chemical,  run='default',
                 chemdb='chemicals_default.txt',
                 seasonalparfile='seasonal_parameters_default.txt',
                 constantparfile='const_parameters_default.txt',
                 compartmentfile='compartments_default.txt',
                 flowdir='default',
                 processfile='processes_default.txt',
                 controlfile='control_default.txt',
                 use_correction = False,       #FY
                 track_flows = False,
                 emissionfile='emissions_default.txt'):    #FY
        ## write model specification to file
        self.track_flows = track_flows
        self.use_correction = use_correction
        self.output_summary(chemical, run, chemdb, seasonalparfile, constantparfile,\
                            compartmentfile, flowdir, processfile, controlfile)       # HW  
        ## Initialize model run
        self.run=run
        
        ## Inintialize Variables for linear flux approximations SSchenker
        self.flux_res=[]
        self.flux_res_error = []
        self.flux_key = {}
        
        #print('\nInitializing model run %s\n') % self.run       
        self.chemical=chemical
        ## read chemical
        fn=os.path.join('Chemicals', chemdb)
        self.update_chemical(chemical, fn)
        ## read environnmental parameters
        fn1=os.path.join('Environment', constantparfile)
        fn2=os.path.join('Environment', seasonalparfile)
        self.update_env_params(fn1,fn2)
        ## read compartment descriptions
        fn=os.path.join('Environment', compartmentfile)
        self.update_compartments(fn)
        ## read flows of transport media
        directory=os.path.join('Flows', flowdir)
        self.update_flows(directory)
        ## read processlist
        fn=os.path.join('Processes', processfile)
        self.update_processlist(fn)
        ## read controlfile
        fn=os.path.join('Control',controlfile)
        self.update_control(fn)
        ## read emissionfile
        self.emissionfile=emissionfile     #FY
        ## update model-characteristic class attributes
        ## nocells, nots, nocomps, matdim
        self.update_char_numbers()

        ## construct time-averaged parameter array and flow-matrices
        ## parmean and flowdictmean
        #self.average_parameters()
        #self.average_flows()
        ## Calculate volumes of compartments and sub-compartments [self.vdict]
        self.update_volumes() # 
        ## update self.killidx, the vector of indices in 0:matdim indicating
        ## zero-volume compartments
        self.mkkillidx_fromV()
        ## Calculate chemical properties in all cells at all timesteps
        self.update_chempardict() #self.chempardict
        ## Calculate z-values [self.zdict]
        self.update_zvalues() 
        ## calculate D-values for intra-cell processes [self.Dproc]
        self.update_dvalues_process()
        ## calculate D-values for flow in transport media [self.Dflow]
        self.update_dvalues_flows()
        ## calculate D-values with correction Q for flow in transport media [self.DflowQ]     Fangyuan 28.10.2016
        self.update_dvalues_flows_Q()
        self.update_dvalues_mix_Q()
        ## make 1/ZV available as list of diagonal sparse
        ## csr matrices [self.ZVinvlist]
        self.update_zvinvlist()
        # Construct the system matrices of D-values,
        # in sparse "coo" format, one for each timestep.
        # Compress matrices (delete 0 volume compartments)
        # and remember positions of deleted rows/columns
        # in self.killidx
        self.update_bigDlist()
    #################### end init ##############################################
        
    ## functions that read input files and construct or update the
    ## class attributes
    ## chemdict, par, compdict, flowdict, proclist, controldict
    def update_chemical(self,chemno,chemdb):
        """reads the chemical to be modelled from the database of chemicals and updates Model.chemdict.      
        :param chemno: ID of the chemical in chemdb.
        :param chemdb: path to the database of chemicals."""
        self.chemdict = readChemicals([chemno], chemdb)[chemno]
    def update_env_params(self,cpfile,spfile):
        self.par = constructEnvironment(cpfile,spfile)
        self.par = modparams(self)
    def update_compartments(self, compfile):
        self.compdict = readCompartments(compfile)
    def update_flows(self, flowdir, mean=False):
        self.flowdict = readFlows(self.compdict.keys(),flowdir)
    def update_processlist(self, processfile):
        self.proclist = readProcesses(self.compdict, processfile)
    def update_control(self, controlfile):
        self.controldict = readControl(controlfile)
    def update_fluxkey(self): # SSchenker
         ''' Constructs the flux dictionary for tracking flux outputs'''
         lat = int(math.sqrt(self.nocells/2))
         lon = lat*2
         use_correction=self.use_correction
         #self.flux_key = mkFluxkey(self.Dflow, self.Dproc, lat=lat, lon=lon)
         self.flux_key = mkFluxkey(self, self.Dflow, self.Dproc, lat=lat, lon=lon, \
                                   use_correction=use_correction)   #Fangyuan 2rd Dec., 2016
    ############################################################################

    ## update of model-characteristic class attributes
    def update_char_numbers(self):
        self.nocells = self.par.shape[0] # number of cells
        self.nots = self.par.shape[1] # number of timesteps in a period
        self.nocomp=len(self.compdict) # number of compartments
        self.matdim=self.nocells*self.nocomp # system-matrix dimension
    ############################################################################

    ## functions that operate on input parameters to
    ## create or update dictionaries:
    ## vdict:{<compartment>:{<subcompartment>:array((nocells, notimesteps))}}
    ## chempardict:{<compartment>:array((nocells,notimesteps),
    ##                                  dtype=[k_reac,Kaw,Koa,Kow])
    ## zdict:{<compartment>:{<subcompartment>:array((nocells, notimesteps))}}
    def update_volumes(self, mean=False):
        if not mean:
            self.vdict = Volumes(self.par,self.compdict).V            
        else:
            self.vdict = Volumes(self.parmean,self.compdict).V        
            
    def update_chempardict(self, mean=False):
        if not mean:
            self.chempardict=tempcorrect(self.par, self.compdict, self.chemdict)
        else:
            self.chempardict=tempcorrect(self.parmean, self.compdict,\
                                         self.chemdict)

    def update_zvalues(self, mean=False):
        if not mean:
            self.zdict = Zvalues(self.par, self.compdict, self.chempardict).Z
        else:
            self.zdict=Zvalues(self.parmean, self.compdict, self.chempardict).Z
    ###########################################################################
    
    ## Calculate D-values 
    def update_dvalues_process(self):
        self.procobj=process(self)
        self.Dproc=self.procobj.getD()
        
    def update_dvalues_flows(self):
        self.Dflow=mkflowD(self)
            
    def update_dvalues_flows_Q(self):                  #FY
        if self.use_correction==True:
            self.DflowQ=mkflowDQ(self)

    def update_dvalues_mix_Q(self):                  #FY
        if self.use_correction==True:
            self.DmixQ=mkDmixQ(self)    
    ############################################################################

    ## Construction of system matrices
    def update_zvinvlist(self):
        '''Make list of 1/ZV diagonal matrices (matdim,matdim)
        in sparse csr-format'''
        print('Making 1/ZV - matrices')
        self.ZVinvlist=mkZVinv(self)
       # print('1/ZV - matrices DONE') #Added by RKG, 2014/06/04

    def mkkillidx_fromV(self):
        #Calculating system-matrix indices indicating
        #zero-volume compartments; setting self.killidx
        killidx=[]
        for c in self.compdict.keys():
            zerovols=where(sum(self.vdict[c]['bulk'], axis=1)==0)[0]
            killidx.extend(tocell(zerovols+1,c,self))
        self.killidx=sort(array(killidx))

    def mkkillidx_fromD(self):
        #Calculating 0-volume compartment indices from Dvalue-matrices.
        #Used for consistency check.
        zerovols=[]
        for m in self.bigDlist:
            zerovols.append(findemptycells(m))
        if not all([x==zerovols[0] for x in zerovols]):
            sys.exit("update_bigDlist: detected inconsistent zero-volumes \
            between timesteps. Aborting !")
        return(array(zerovols[0]))

    def _DtoK(self):
        #Create mass rate-constant matrices from D-value matrices
        for m in arange(0,len(self.bigDlist)):
            self.bigDlist[m]= self.bigDlist[m]*(self.ZVinvlist[m])
    
    def update_bigDlist(self):
        '''Make system matrices ready for solving'''
        track_flows=self.track_flows
        use_correction=self.use_correction
        self.bigDlist=mkbigD(self, track_flows=track_flows, use_correction=use_correction)
        if self.controldict['dumpDmatrices'] in ['1','y','Y','yes','Yes','YES']:
            fn=os.path.join('Output', self.run, 'matrixdump.cpk')
            if not os.path.exists(os.path.dirname(fn)):
                os.mkdir(os.path.dirname(fn))
            else:
                if os.path.exists(fn):
                    print('Attention: overwriting %s\n') % (fn)
            of=open(fn,'w')
            cPickle.dump(self.bigDlist,of,protocol=2)
            of.close()
            print('wrote D-value - matrices to %s\n') % (fn)
        if self.controldict['dumpDmatricesTxt'] in ['1','y','Y','yes','Yes','YES']:   # txt-output added by HW, 15.01.2011
            fn2=os.path.join('Output', self.run, 'matrixdump.txt')
            fn3=os.path.join('Output', self.run, 'matrixtime.txt')   #FY
            if not os.path.exists(os.path.dirname(fn2)):
                os.mkdir(os.path.dirname(fn2))
            else:
                if os.path.exists(fn2):
                    print('Attention: overwriting %s\n') % (fn2)
            if (not track_flows):
                write_bigDlist_txt(self.bigDlist, fn2, fn3,self) #FY                
                print('wrote D-value - matrices to %s\n') % (fn2)
                print('wrote D-value - matrices to %s\n') % (fn3)
                # end of modification

        
        ## consistency check for killidx

        print('Consistency check of zero-volume compartments: '),
        kidx=self.mkkillidx_fromD()
        if track_flows:
           killidx_double = concatenate(\
                   (self.killidx, add(self.killidx, self.matdim)),\
                   axis=0)
        else: 
            killidx_double = self.killidx
           
        if not all(kidx == killidx_double):
            sys.exit('Detected inconsistencies between\
            0-volumes and 0 rows/columns in D-matrices')
        print(' passed\n')
        
        ## compress D-value list
        for m in arange(0,len(self.bigDlist)):
            self.bigDlist[m]=compressmat(self.bigDlist[m], killidx_double)
                   
        ## transform bigDlist to mass rate constants (divide by ZV)
        ## and save in csr - sparse format
        if track_flows:
            for ii in range(len(self.ZVinvlist)):
                self.ZVinvlist[ii] =\
                    sp.bmat([[self.ZVinvlist[ii], None],\
                            [None, self.ZVinvlist[ii]]])

                
        self._DtoK()
    ############################################################################

    ## Functions reading input that are not called in __init__
    def update_emissions(self,emfile):
        ''' Constructs emission object from emissions-input-file'''
        fn=os.path.join('Emissions', emfile)
        self.emission = Emission(fn)
        self.output_summary2(self.run, emfile)  # HW
        
    def update_solver(self, stepfile='stepping_default.txt',
                      solvparamfile='solvparams_default.txt',
                      initfile='initial_default.txt'):
        '''
        Updates or initializes the solver.

        :param string stepfile: File describing the mapping
          between timesteps and actual time
        :param string solvparamfile: File containing the parameters
          for the solver
        :param string initfile: File containing the initial
          contamination of the environments.
        '''
        stepfile=os.path.join('Solver', stepfile)
        solvparamfile=os.path.join('Solver', solvparamfile)
        initfile=os.path.join('Solver', initfile)
        self.solver=Solver(self,stepfile,solvparamfile,initfile)
        self.output_summary3(self.run, solvparamfile)  # HW


    ## Functions that solve the system
    def solve_ss_original(self): # This is the original solve_ss
        '''
        Averages the all mass rate matrices for each season (month) to
        get a time invariant system. Solves for the steady-state.
        '''
        ## average D-value matrices
        self.Davg = mean(self.bigDlist)
        try:
            q=self.emission.get_emission(0,self,type='array')
        except AttributeError:
            print("update_emissions(<filename>) has to be called before\n"\
                  +"solving for steady-state. Did nothing!\n")
            return()
        ss_res =dot(-inv(self.Davg.toarray()), q)
        self.ss_res=expandvec(ss_res, self.killidx)

    def solve_ss(self): # This is the modified solve_ss, RKG, 24.07.2014
        '''
        Averages the all mass rate matrices for each season (month) to
        get a time invariant system. Solves for the steady-state.
        '''
        ## average D-value matrices
        #self.Davg = mean(self.bigDlist) # gives an error, RKG, 24.07.2014
        # RKG's code starts here, 24.07.2014
        bigDlist_temp = []
        for i in range(self.nots):
            bigDlist_temp.append(self.bigDlist[i].toarray())
        #bigDarray = array(bigDlist_temp)
        #Davg_temp = mean(bigDarray, axis = 0)
        Davg_temp = mean(bigDlist_temp, axis = 0)
        self.Davg = sp.coo_matrix(matrix(Davg_temp))
        # RKG's code ends here, 24.07.2014
        try:
            q=self.emission.get_emission(0,self,type='array')
        except AttributeError:
            print("update_emissions(<filename>) has to be called before\n"\
                  +"solving for steady-state. Did nothing!\n")
            return()
        ss_res =dot(-inv(self.Davg.toarray()), q) # self.Davg is already an array now, RKG, 24.07.2014
        #ss_res =dot(-inv(self.Davg), q) # line added by RKG, 27.04.2014
        self.ss_res=expandvec(ss_res, self.killidx)
    
    def solve_dyn(self, no_periods, with_deriv = False, use_odespy=False):
        '''
        Solves the dynamic model for *no_periods* (years).
        '''
        if with_deriv or self.track_flows:
            self.update_fluxkey()
        if with_deriv:
            if use_odespy:
                self.dyn_res, self.flux_res, self.flux_res_error = \
                     self.solver.solve_ode_odespy(no_periods,
                                                  with_deriv=with_deriv)
                
            else:
                self.dyn_res, self.flux_res, self.flux_res_error = \
                     self.solver.solve_ode(no_periods,
                                           with_deriv=with_deriv)
        else:
            if use_odespy:
                self.dyn_res, self.flux_res = \
                     self.solver.solve_ode_odespy(no_periods,
                                                  with_deriv=with_deriv,
                                                  track_flow=self.track_flows)
            else:
                self.dyn_res, self.flux_res = \
                     self.solver.solve_ode(no_periods,
                                           with_deriv=with_deriv,
                                           track_flow=self.track_flows)

        reslist=[]
        #print 'solve_dyn: self.dyn_res:', type(self.dyn_res),self.dyn_res.shape # RKG, 03.07.2014
        for resvec in self.dyn_res.T:
            reslist.append(expandvec(resvec, self.killidx))
        self.dyn_res=array(reslist).T
        #print 'solve_dyn: self.dyn_res:', type(self.dyn_res),self.dyn_res.shape # RKG, 03.07.2014
        
        
        for i in range ( len(self.flux_res)): #SSchenker Linear Flux Approx. 
            #print "data processing flux from season %i" % i
            print("data processing flux from season %i" % i) # modification to run with python 3, RKG, 10.06.2014 RKG 
            self.flux_res[i] = self.flux_res[i].tolil()
            
            
    ############################################################################

    ## output routines
    def output_summary(self, nchemical, nrun, nchemdb, nseasonalparfile, nconstantparfile,\
                       ncompartmentfile, nflowdir, nprocessfile, ncontrolfile):
        fn=os.path.join('Output', nrun, 'summary.txt')
        write_output_summary(self, fn, nchemical, nrun, nchemdb, nseasonalparfile, nconstantparfile,\
                       ncompartmentfile, nflowdir, nprocessfile, ncontrolfile) 
        
    def output_summary2(self, nrun, nemfile):
        fn=os.path.join('Output', nrun, 'summary.txt')
        write_output_summary2(self, fn, nemfile) 
        
    def output_summary3(self, nrun, nsolvfile):
        fn=os.path.join('Output', nrun, 'summary.txt')
        write_output_summary3(self, fn, nsolvfile) 
    
    def output_ss(self, filename='ss', units=['mol'], netcdf=False, cpk=True):  # parameter cpk added by HW
        """
        :param string filename: Prefix of the output file(s). The output files
                                will be written to Output/<run>/<filename>_x.
        :param list units: A list of strings to indicate the unit(s) of the
         output. Possible list elements are 'mol', 'mol_per_m3', 'kg',
         'kg_per_m3' and 'Pa'.
        :param bool netcdf: Indicates whether netCDF-files are written or not.
        """
        ## filename without extension (automatically generated from filetype)
        fn=os.path.join('Output',self.run,filename)
        write_output_ss(self,fn,units,netcdf,cpk)

    def output_ss_txt(self, filename='ss'):  # Function added by RKG, 04.08.2014
        """
        :param string filename: Prefix of the output file(s). The output files
                                will be written to Output/<run>/<filename>_x.
        Currently, the unit of the output is 'mol'.
        """
        ## filename without extension (automatically generated from filetype)
        fn=os.path.join('Output',self.run,filename)
        write_output_ss_txt(self,fn)

    def output_dyn(self,filename='dyn',units=['mol'], netcdf=False, cpk=True):   # parameter cpk added by HW
        """
        :param string filename: Prefix of the output file(s). The output files
                                will be written to Output/<run>/<filename>_x.
        :param list units: A list of strings to indicate the unit(s) of the
         output. Possible list elements are 'mol', 'mol_per_m3', 'kg',
         'kg_per_m3' and 'Pa'.
        :param bool netcdf: Indicates whether netCDF-files are written or not.
        """
        ## filename without extension (automatically generated from filetype)
        fn=os.path.join('Output',self.run,filename)
        write_output_dyn(self,fn,units,netcdf,cpk)
        
    def output_se(self,
                  filename='sec_em',
                  units=['mol'],
                  netcdf=False,
                  cpk=True,
                  scenario="air"):   # parameter cpk added by HW
        """
        :param string filename: Prefix of the output file(s). The output files
                                will be written to Output/<run>/<filename>_x.
        :param list units: A list of strings to indicate the unit(s) of the
         output. Possible list elements are 'mol', 'mol_per_m3', 'kg',
         'kg_per_m3' and 'Pa'.
        :param bool netcdf: Indicates whether netCDF-files are written or not.
        """
        ## filename without extension (automatically generated from filetype)
        fn=os.path.join('Output',self.run,filename)
        write_output_se(self,fn,units,netcdf,cpk, scenario=scenario)       
        
    def output_fluxes(self,filename='fluxes',units=['mol'], netcdf=False):
        """
        SSchenker: (Linerar Approx. Fluxes)
        :param string filename: Prefix of the output file(s). The output files
                                will be written to Output/<run>/<filename>_x.
        :param list units: A list of strings to indicate the unit(s) of the
         output. Possible list elements are 'mol', 'mol_per_m3', 'kg',
         'kg_per_m3' and 'Pa'.
        :param bool netcdf: Indicates whether netCDF-files are written or not.
        """
        ## filename without extension (automatically generated from filetype)
        lat = int(math.sqrt(self.nocells/2))
        lon = lat*2
        fn=os.path.join('Output',self.run,filename)
        write_output_dflux(self,fn,
                           units=['mol'],
                               netcdf=netcdf,
                           nlat=lat,
                           nlon=lon)        
        
    def output_end_txt(self, filename='endstate.txt', tstep=13):                         # added by HW
        fn=os.path.join('Output',self.run,filename)
        write_output_end_txt(self,fn,tstep)
                   
    def output_dyn_tmp(self, filename='initial_tmp.txt', tstep=13):                        # added by HW
        fn=os.path.join('Solver', filename)
        write_output_dyn_tmp(self,fn,tstep)
        
    ############################################################################

    ## accessors
    def get_avgMatrix(self):
        """
        Returns the averages system matrix of mass rate constants used
        for the steady-state solution.
        """
        print("Constructing full steady-state-matrix ...")
        try:
            avgmatrix=expandmatrix(self.Davg, self.killidx).toarray()
        except AttributeError:
            print("solve_ss() has to be called before acessing the averaged"\
                  +" system matrix.")
            return()
        print("OK")
        return(avgmatrix)

    def get_ssEmission(self):
        """
        Returns the full emission vector used for steady-state
        calculation in Compressed Sparse Column format.
        """
        try:
            ssem=self.emission.get_emission(0,self,'array')
        except AttributeError:
            print("update_emissions(<filename>) has to be called before"\
                  +" acessing the emission vector. Did nothing!")
            return()
        print("Expanding emission vector ...")
        ssem=expandvec(ssem, self.killidx)
        print("OK")
        return(ssem)
    
    ## construct or update time-averages of parameters and flow-matrices
    def average_parameters(self):
        self.parmean=empty((self.par.shape[0]), dtype=self.par.dtype)
        for field in self.par.dtype.names:
            self.parmean[field]=mean(self.par[field], axis=1)
    def average_flows(self):
        self.flowdictmean={}
        for item in self.flowdict.items():
            meanmat=c_[item[1][:,0:2], mean(item[1][:,2:], axis=1)]
            self.flowdictmean[item[0]]=meanmat

############################# END OF BETRS.py ##################################

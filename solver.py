################################################################################
'''Module "solver.py", July 2010, is part of
BETR-Research by Harald von Waldow <hvwaldow@chem.ethz.ch>, which is
based on BETR-Global by Matt MacLeod <matthew.macleod@chem.ethz.ch>

This module contains the solver for BETR-Research'''
################################################################################
from numpy import *
import inspect
import pdb
import sys
import re
import math
import os.path
import time
from scipy.sparse import *
try:
    from scipy.integrate import ode
except  ImportError:
    from scipy.integrate.ode import ode
try:
    import odespy
except  ImportError:
    pass
from warnings import warn

from scipy.linalg import *
from helpers import *
import pickle

class Solver:
    def __init__(self, model, stepfile, solvparamfile, initfile):
        self.model=model
        self.y0=self._readinitial(model,initfile)
        self.solvparms=self._readsolvparams(solvparamfile)
        self.stepdef, self.period, self.steps_per_period=\
                      self._readsteps(stepfile)
        
    def _readinitial(self,model,filename):
        y0=zeros((model.matdim))
        f=open(filename,'r')
        while True:
            line=f.readline()
            if line == '': break
            if line[0] =='#' or line == '\n': continue
            line=line.split()
            cell=tocell(int(line[0]), int(line[1]),model)
            y0[cell]=float(line[2])
        f.close()
        y0=compressvec(y0,model.killidx)
        return(y0)
    
    def _readsolvparams(self,filename):
        paramdict={}
        f=open(filename,'r')
        while True:
            line=f.readline()
            if line == '': break
            if line[0] =='#' or line == '\n': continue
            line=line.split()
            paramdict[line[0]]=line[1]
        f.close()
        for p in paramdict:
            try:
                paramdict[p]=int(paramdict[p])
            except ValueError:
                try:
                     paramdict[p]=float(paramdict[p])
                except ValueError: continue
        return(paramdict)

    def _readsteps(self, stepfile):
        f=open(stepfile, 'r')
        lines=f.readlines()
        f.close
        for lin in lines:
            if lin[0] == '#' or lin[0] == '\n':
                continue
            break
        stepdef = array([float(x) for x in lin.split()])
        period = sum(stepdef)
        steps_per_period =len(stepdef)
        return (stepdef, period, steps_per_period)
        
    def get_stepidx(self,t):
        """ Get the step that belongs to t. Steps start with 0
        returns (absolute stepidx, step_in_period)"""
        t=float(t)
        if (t == 0):
            return (0, 0)
        # periodidx = floor(t / self.period)
        t_in_period = mod(t,self.period)
        for cs in enumerate(cumsum(self.stepdef)):
            if t_in_period <= cs[1]:
                step_period = cs[0]
                break
        stepidx = floor(t/self.period)*self.steps_per_period+step_period
        return (int(stepidx), int(step_period))
    

    def solve_ode(self, no_periods, with_deriv=False, exp_t_steps=False, track_flow=False):
        """Set basic parameters"""
        ts = self.model.nots*no_periods
        
        res = self.y0
    
        if track_flow:
            self.y0 =concatenate((self.y0 ,zeros(shape(self.y0))),axis=0)
            
            
        
        flux_res = []
        flux_res_error = []
        
        """Set parameters for ode integrtation"""
        ig = ode(f,jac)
        
        ig.set_integrator(self.solvparms['integname'],\
                          atol = self.solvparms['atol'],\
                          rtol = self.solvparms['rtol'],\
                          method = self.solvparms['method'],\
                          with_jacobian = self.solvparms['with_jacobian'],\
                          order = self.solvparms['order'],\
                          nsteps=self.solvparms['nsteps']\
                          )

        ig.set_initial_value(self.y0,0)
        """ Start timer"""
        ts0=time.time()
        
        for m in range(0, ts):
            print('Integrating month # %s') % (m)
            #ig.set_initial_value(ig.y, ig.t)
            
            current_season = mod(m, self.steps_per_period)
            current_matrix_csr = self.model.bigDlist[current_season]
            current_emissions = self.model.emission.get_emission(m,self.model)
            #current_Q = self.model.DflowQ[current_season]  #Fangyuan  22rd Nov. 2016
            if track_flow:
                 current_emissions =\
                     csc_matrix(vstack((current_emissions.todense(),\
                     zeros(shape(current_emissions)))))
#                fx=open("CE.pyc", "w")
#                pickle.dump(current_emissions, fx)
#                fx.close()
            
            ig = ig.set_f_params(current_matrix_csr, current_emissions)
            current_jac = current_matrix_csr.toarray()
            ig = ig.set_jac_params(current_jac)
            no_sub_steps = int(ceil(self.stepdef[current_season]
                                    /self.solvparms['maxsubsteplength']))
                                 
            
            if with_deriv and not(exp_t_steps) and no_sub_steps<146:
                no_sub_steps=146
                substeplength = self.stepdef[current_season] / no_sub_steps
                substeplist = ones(no_sub_steps) * substeplength

                Diff_upper = (dot(self.model.bigDlist[current_season],
                                                    diags(self.y0, 0))
                                                    * substeplength)
                print('\tSeason %d\t %d substeps of length %.5f hours')\
                            % (current_season,no_sub_steps, substeplength)
                                   
            elif with_deriv and exp_t_steps:
                substeplist= (10.**(1./10.))**arange(1,11)
                substeplist=\
                    (substeplist*self.stepdef[current_season]) \
                    / sum(substeplist)
                
                Diff_upper = (dot(self.model.bigDlist[current_season],
                                                    diags(self.y0, 0))
                                                    * substeplist[0])
                no_sub_steps=len(substeplist)                                                   
                print('\tSeason %d\t %d substeps of length %.5f to %.5f hours')\
                            % (
                               current_season,
                               no_sub_steps,
                               substeplist[0],
                               substeplist[-1],
                               )                                                    
                
            else:    
                substeplength = self.stepdef[current_season] / no_sub_steps

                print('\tSeason %d\t %d substeps of length %.5f hours')\
                            % (current_season,no_sub_steps, substeplength)
            te1=time.time()
            ts1=time.time()
            print('\tIntegrating %i substeps: ' % no_sub_steps),
            for i in range(0, no_sub_steps):

                ig.integrate(ig.t + substeplength)
                if with_deriv:
                    if i == 0:                        
                        '''
                        In the moment the full flux matrix is calculated. 
                        It would be absolutely sufficient to calculate only the
                        diagonal elements, which already contain all of the
                        information 
                        #SSchenker
                        '''
                        Diff_lower = (dot(self.model.bigDlist[current_season],
                                                    diags(ig.y, 0))
                                                    * substeplist[i])
                    else:
                        Diff_lower = (Diff_lower + 
                                     (dot(self.model.bigDlist[current_season],
                                     diags(ig.y, 0))
                                     * substeplist[i])
                                     )
                    
                    if i != (no_sub_steps-1):
                        Diff_upper = (Diff_upper + 
                                     (dot(self.model.bigDlist[current_season],
                                     diags(ig.y, 0))
                                     * substeplist[i+1])
                                     )
                
                if (not ig.successful()):
                    sys.exit("Integration failed; aborting !")
                te1=time.time()
                if no_sub_steps > 1:
                    if i % (no_sub_steps/20) == 0: print('.'),
            
            print('\n done with month %d  [%.3f s]')  % (m, (te1-ts1))
            
            
            if track_flow:
                y_new = split(ig.y, 2)
                res=vstack((res, y_new[0]))
                flux_res.append(lil_matrix(y_new[1].T * eye(len(y_new[1]), len(y_new[1])) ))
                
                
                y_new = concatenate((y_new[0] ,zeros(shape(y_new[0]))),axis=0)
                ig = ode(f,jac)
        
                ig.set_integrator(self.solvparms['integname'],\
                          atol = self.solvparms['atol'],\
                          rtol = self.solvparms['rtol'],\
                          method = self.solvparms['method'],\
                          with_jacobian = self.solvparms['with_jacobian'],\
                          order = self.solvparms['order'],\
                          nsteps=self.solvparms['nsteps']\
                          )

                ig.set_initial_value(y_new,0)
                
            elif with_deriv:
                res=vstack((res, ig.y))
                flux_res.append(0.5 * Diff_lower + 0.5 * Diff_upper)
                flux_res_error.append((Diff_upper-Diff_lower)/2)
            else:
                #ig.set_initial_value(ig.y*current_Q.toarray()[0], ig.t)   #Fangyuan 22rd Nov. 2016
                #res=vstack((res, ig.y*current_Q.toarray()[0]))      #Fangyuan 22rd Nov. 2016
                res=vstack((res, ig.y))
            """Stop timer"""
            te0=time.time()
        print "time elapsed: %f minutes" % ((te0-ts0)/60.0)
        res=res.T
        if track_flow:
            return res, flux_res
        elif with_deriv:
            return res, flux_res, flux_res_error
        else:
            return res, []
            
    def solve_ode_odespy(self, no_periods, with_deriv=False, exp_t_steps=False, track_flow=False):
        """
        Solve the Sytem of ODEs with the ODESPY Runge-Kutta-Fehlberg 45
        needs odespy to function properly
        """
        if "odespy" not in sys.modules:
            warn("""
ODESPY can not be imported: proceeding with the interal ODE solver:
please install ODESPY for increased performance
or set use_odespy=false
                """)  
            return
        
        
        
        '''Check for unsupported options '''
        if with_deriv:
            raise NotImplementedError

        '''Initialize basic variables'''            
        ts = self.model.nots*no_periods
        res = self.y0
        current_y = self.y0
        
        if track_flow:
            current_y =concatenate((self.y0 ,zeros(shape(self.y0))),axis=0)
        flux_res = []

        """ Start timer"""
        ts0=time.time()
        
        '''Cycle through months m'''
        for m in range(0, ts):
            print('Integrating month # %s') % (m)
            current_season = mod(m, self.steps_per_period)
            current_matrix_csr = self.model.bigDlist[current_season]
            current_emissions = self.model.emission.get_emission(m,self.model)

            if track_flow:
                 current_emissions =\
                     csc_matrix(vstack((current_emissions.todense(),\
                     zeros(shape(current_emissions)))))

            """Set parameters for ode integrtation"""
            t=0
            mb = MassBalance(current_matrix_csr, current_emissions)
            solver = odespy.RKF45(mb.f,  
                                        atol = self.solvparms['atol'],\
                                        rtol = self.solvparms['rtol'],\
                                        )
            solver.set_initial_condition(current_y)
            no_sub_steps = 1
            substeplength = self.stepdef[current_season] 
            print('\tSeason %d\t %d substeps of length %.5f hours')\
                            % (current_season,no_sub_steps, substeplength)
            te1=time.time()
            ts1=time.time()
            print('\tIntegrating %i substeps: ' % no_sub_steps),
            
            '''NOTE: Check if loopng over substeps cam be avoided (should always be 1 step) '''
            for i in range(0, no_sub_steps):
                out_y , t = solver.solve([t,substeplength+t])     #t =[0,730], substeplength=730
                out_y=out_y[1]
                t=t[1]    
                te1=time.time()
            
            print('\n done with month %d  [%.3f s]')  % (m, (te1-ts1))
            
            
            if track_flow:
                '''Upper part of y are c(X), lower part are integrated mass fluxes'''
                current_y = split(out_y, 2)
                res=vstack((res, current_y[0]))
                flux_res.append(lil_matrix(current_y[1].T * eye(len(current_y[1]), len(current_y[1])) ))
                current_y = concatenate((current_y[0] ,zeros(shape(current_y[0]))),axis=0)
            else:
                '''Only "upper" part of y exists  '''
                current_y=out_y
                res=vstack((res, out_y))
               #out_y=out_y*current_Q.toarray()[0]   #Fangyuan  22rd Nov. 2016
                #current_y=out_y            #Fangyuan  22rd Nov. 2016
                #res=vstack((res, out_y))   #Fangyuan  22rd Nov. 2016
            """Stop timer"""
            te0=time.time()

        print "time elapsed: %f minutes" % ((te0-ts0)/60.0)
        res=res.T
        
        if track_flow:
            return res, flux_res
        else:
            return res, []            
            
# ##### end class Solver #######################################################

def f(t, y, A, em):
    y = csr_matrix(matrix(y)).T
    dydt = A * y + em
    dydt = dydt.toarray().ravel()
    return(dydt)
    
def fort_f(A,em):
    f_f77_str = """
    subroutine f_f77(neq, t, u, udot)
    Cf2py intent(hide) neq
    Cf2py intent(out) udot
    integer neq
    double precision t, u, udot
    dimension u(neq), udot(neq)
    udot(1) = %.3f*u(1)*(1 - u(1)/%.1f)
    return
    end
    """ % (a, R)
    return f_f77_str
    

def jac(t, y, jac):
    return(jac)
    
################################################################################
class MassBalance:
    def __init__(self, A, em):
        self.em = em
        self.A = A
    
    def f(self, y, t):
        y = csr_matrix(matrix(y)).T
        dydt = self.A * y + self.em
        dydt = dydt.toarray().ravel()
        return(dydt)


#    def u_exact(self, t):
#        a, R, A = self.a, self.R, self.A  # short form
#        return R*A*exp(a*t)/(R + A*(exp(a*t) - 1))

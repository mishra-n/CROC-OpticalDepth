import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.stats as stats
import scipy.interpolate as interp
import scipy.integrate as integrate
import scipy as sp
import re
from matplotlib import rc
import scipy.ndimage as ndimage
from matplotlib import cm
import os
import time as time
import camb
from numba import jit

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit

np.random.seed(2021)
rc('font',**{'size':'20','family':'serif','serif':['CMU serif']})
rc('mathtext', **{'fontset':'cm'})
rc('text', usetex=True)
rc('legend',**{'fontsize':'13'})

CROC_PARAMS ={  'H0':      68.14,
                'OmegaM':  0.3036,
                'OmegaD':  0.2557,
                'OmegaB':  0.0479,
                'OmegaL':  0.6964,
                'OmegaK':  0,
                'DeltaDC': 0.0722008}        

class Universe:
    def __init__(self, h=0.7, omega_b=0.046, omega_m=0.286, sigma8=0.82, n_s=0.96, w=-1, tau=0.07, param_dict=None):
        self.c = 3e10 #cm/s 
        self.sig_T = 6.65e-25 #cm^2

        self.h = h
        self.omega_b = omega_b
        self.omega_m = omega_m
        
        if param_dict is not None:

            self.h = param_dict['H0']/100
            self.omega_b = param_dict['OmegaB']
            self.omega_m = param_dict['OmegaM']


        self.sigma8 = sigma8
        self.n_s = n_s
        self.w = w
        self.tau = tau
        
        #derived parameters
        self.H_0 = h*100 #km/s/Mpc
        self.omega_cdm = self.omega_m - self.omega_b
        
        #little h parameters
        self.ombh2 = self.omega_b * self.h**2
        self.omch2 = self.omega_cdm * self.h**2
        self.ommh2 = self.omega_m * self.h**2
        
    def runCAMB(self):
        
        #Set up a new set of parameters for CAMB
        #pars = camb.CAMBparams()
        pars = camb.set_params(H0=self.H_0, ombh2=self.ombh2, omch2=self.omch2, mnu=0.06, omk=0, tau=self.tau, As=2e-9, ns=self.n_s, r=0, lmax=2500, lens_potential_accuracy=0, WantTransfer=True)
    
        self.CAMBparams = pars

#       #This function sets up CosmoMC-like settings, with one massive neutrino and helium set using BBN consistency
#        pars.set_cosmology(H0=self.H_0, ombh2=self.ombh2, omch2=self.omch2, mnu=0.06, omk=0, tau=self.tau)
#        pars.InitPower.set_params(As=2e-9, ns=self.n_s, r=0)
#        pars.set_for_lmax(2500, lens_potential_accuracy=0);

        #calculate results for these parameters
        results = camb.get_results(pars)
        
        self.CAMB_results = results
    
    def getCAMBPowerSpectrum(self, var1=2, var2=2, hubble_units=True, k_hunit=True, have_power_spectra=True, nonlinear=False, kmax=20, z_range=np.array([0])):
        
        self.CAMBparams.set_matter_power(redshifts=z_range, kmax=kmax)
                
        k,z, P_kz = self.CAMB_results.get_linear_matter_power_spectrum(var1=var1, var2=var2, hubble_units=hubble_units, k_hunit=k_hunit, have_power_spectra=have_power_spectra, params=self.CAMBparams, nonlinear=nonlinear)
        
#        self.k = k
#        self.z = z
#        self.P_kz = P_kz
        
        return k,z,P_kz
        
        
        
        
        
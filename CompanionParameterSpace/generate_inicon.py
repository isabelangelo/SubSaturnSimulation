import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import interpolate
from Kepler1656_companion_parameter_space import *

# ===========================================================
"""
Functions for generating initial conditions of KOI0367 for OSPE

Initial conditions are generated by initial_conditions()
into the following format: 
    
    array([m1,m2,m3,R1,R2,spin1P,spin2P,beta,beta2,gamma,gamma2,a1,a2,e1,e2,g1,g2,i,age])

these initial conditions can be fed into generate_triple.in to generate the initial 
condition files that OSPE uses
"""

a_subsat = 0.197 # separation of subsaturn, au
m_subsat = 1.4585e-4 # 48.6 Mearth (Brady et al. 2018)

# ===========================================================


def draw_sample_companion():
    """
    Function to draw sma and mass of outer companion 
    from parameter space allowed by observations
    :return: a2 in au, m3 in solar mass (tuple)
    """
    # define functions that set upper limits, you want to do this outside of the loop!
    f_RV = interpolate.interp1d(RV_sma, RV_mass)
    f_AO = interpolate.interp1d(AO_sma, AO_mass)
    f_GAIA = 1
    
    # define minima
    min_a2 = 10*a_subsat # au
    min_m3 = 0.1*m_subsat # require mass fraction of 10 for Kozai
    
    # compute a2 : min=a_subsaturn, max from Winn (2015) 
    a2 = np.power(10, np.random.uniform(np.log10(min_a2), 4)) # log uniform [0.1,1e5]
    
    # compute m3 : min=0.1m2, max defined by constraints
    if a2 < AO_sma[0]:
        # RV constraint
        max_m3 = f_RV(a2)
    elif a2 > AO_sma[-1]:
        # GAIA constraint
        max_m3 = 1
    else:
        # whichever is lowest between RV and AO
        max_m3 = min(f_RV(a2), f_AO(a2))
        
    m3 = np.power(10, np.random.uniform(np.log10(min_m3), np.log10(max_m3)))
    
    return(a2,m3)

# ===========================================================

def get_stable_system(a1):
    """
    Function to generate companion mass/sep/ecc 
    and test if stable. Stability criteria requires epsilon < 0.1
    where epsilon = (a1/a2)*(e2/(1-e2**2.))
    :param a1: initial separation of inner planet in au (float)
    :return: a2 in au,m3 in solar mass, e2 (tuple)
    """
    system_found = False
    count = 0
    max_count = 500
    epsilon_max = 0.1
    while system_found == False:
        
        # draw mass [Msun], sma [au] of outer planet from allowed parameter space
        a2, m3 = draw_sample_companion()
        
        # draw eccentricity of outer planet
        e2 = np.random.uniform(0,1)
        
        # compute epsilon for stability criteria
        epsilon = (a1/a2)*(e2/(1-e2**2.))
        
        # apply EKL criteria (close+massive)
        if (m3 < 0.1*m_subsat) and (a2 > 30):
            count += 1
            #print('EKL not satisfied ')
            if count == max_count:
                #print('doesn\'t induce EKL')
                return None
            
        # apply stability criterion
        elif epsilon > 0.1:
            count += 1
            #print('stability not satisfied ')
            if count == max_count:
                #print('couldn\'t find stable system')
                return None
            
        # save systems that satisfy criteria
        else:
            system_found = True
            #print('system found ', count)
            return(a2,m3,e2)
        
# ===========================================================

def initial_conditions():
    """
    Function to generate generate the initial conditions 
    to go in the triple.in file
    :return: array of initial condition values (np.array), 
             or None if no stable system is found
    """
    # masses [Msun]
    m1 = 1.03
    m2 = m_subsat

    # radii [Rsun]
    R1 = 1.1
    R2 = 4.594e-2 #5.02 Rearth (Brady et al. 2018)

    # spins [fraction of day]
    spin1P = 10 
    spin2P = 10 

    # betas [deg]
    beta = 0
    beta2 = 0

    # gamma [deg]
    gamma = 45
    gamma2 = 45 
    
    # g values [deg]
    g1 = np.random.uniform(0,360)
    g2 = np.random.uniform(0,360)

    # mutual inclination [deg]
    cosi = np.random.uniform(np.cos(0),np.cos(np.pi))
    i = np.rad2deg(np.arccos(cosi))

    # age to evolve to [Myr]
    age =  6.31*1.2*1000 # evolve to 1.2*system age (Brady et al. 2018)
    
    # inner planet sma [au]
    a1_min = 0.1 # minimum set by inner radius of dust in protoplanetary disk
    a1 = np.power(10, np.random.uniform(np.log10(a1_min), 1)) # log uniform [0.8*a_subsat,10]
    
    # inner planet eccentricity
    e1 = 0.01
    
    # outer planet sma/mass/ecc
    stable_system = get_stable_system(a1)
    if stable_system is not None:
        a2,m3,e2 = stable_system
        return(np.array([m1,m2,m3,R1,R2,spin1P,spin2P,beta,beta2,gamma,gamma2,a1,a2,e1,e2,g1,g2,i,age])) 
    else:
        return None
    
# ===========================================================
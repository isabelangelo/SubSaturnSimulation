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
    f_AO = interpolate.interp1d(AO_sma, AO_mass)
    f_Speckle = interpolate.interp1d(Speckle_sma,Speckle_mass)
    f_GAIA = 1
    
    # define minima
    min_a2 = 1.40 # au
    min_m3 = 0.1*m_subsat # require mass fraction of 10 for Kozai
    
    # compute a2 : 600 days < P < 2100 days from RV data
    a2 = np.random.uniform(1.40,3.25) # uniform [1.40au,3.25au] 
    
    # compute m3 : min=0.1m2, max defined by constraints
    if a2 < AO_sma[0] or a2 > AO_sma[-1]:
    	# maximum set by theory to 2Mjup
    	max_m3 = 0.0019 # Msun
    elif a2 < Speckle_sma[0] or a2 > Speckle_sma[-1]:
    	# AO constraint
    	#max_m3 = f_AO(a2) from graph
        max_m3 = 0.0019 # Msun # from theory
    else:
    	# whichever is lowest between Speckle and AO
    	max_m3 = min(f_Speckle(a2), f_AO(a2))
        
    m3 = np.power(10, np.random.uniform(np.log10(min_m3), np.log10(max_m3)))
    
    return(a2,m3)

# ===========================================================

def get_stable_system(a2,m3,e2):
    """
    Function to generate companion mass/sep/ecc 
    and test if stable. Stability criteria requires epsilon < 0.1
    where epsilon = (a1/a2)*(e2/(1-e2**2.))
    :param a2,m3,e2: sma, mass, ecc of outer planet (au, Msun, '')
    :return: a1 in au
    """
    system_found = False
    count = 0
    max_count = 500
    epsilon_max = 0.1
    while system_found == False:

        # inner planet sma [au]
        a1_min = 0.1 # minimum set by inner radius of dust in protoplanetary disk
        a1_max = epsilon_max*a2*((1-e2**2.)/e2) # set by epsilon<0.1 for stability
        a1 = np.power(10, np.random.uniform(np.log10(a1_min), np.log10(a1_max)))

        # compute epsilon for stability criteria
        epsilon = (a1/a2)*(e2/(1-e2**2.))

        # apply EKL criteria (close+massive)
        if (m3 < 0.1*m_subsat) and (a2 > 30):
            count += 1
            #print('EKL not satisfied ')
            if count == max_count:
                #print('doesn\'t induce EKL')
                return None

        # apply stability criterion (epsilon>0.1 means a1max<a1min)
        elif epsilon > epsilon_max:
            count += 1
            #print('stability not satisfied ')
            if count == max_count:
                #print('couldn\'t find stable system')
                return None

        # assert that a1<a2, reject if not
        elif a2 < a1:
            count += 1
            #print('a2 < a1')
            if count == max_count:
                #print('couldn\'t find stable system')
                return None
            
        # save systems that satisfy criteria
        else:
            if a2 < a1:
                print('a2<a1')
            system_found = True
            #print('system found ', count)
            return(a1)
        
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

    # outer planet sma/mass/ecc
    a2, m3 = draw_sample_companion()
    e2 = np.random.uniform(0,1)

    # inner planet eccentricity
    e1 = 0.01
    
    # outer planet sma/mass/ecc
    stable_system = get_stable_system(a2,m3,e2)
    if stable_system is not None:
        a1 = stable_system
        return(np.array([m1,m2,m3,R1,R2,spin1P,spin2P,beta,beta2,gamma,gamma2,a1,a2,e1,e2,g1,g2,i,age]))
    else:
        return None


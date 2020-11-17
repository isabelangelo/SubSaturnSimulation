
"""
This code converts observational data into constraints
on the parameter space in which a companion to Kepler-1656 
could reside.
"""

from cgs import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from arcsec2au import *
from mag2mass import *

system_distance = 186.9 # pc


# ======== AO data ===========================================
# sep in arcsec, dmag in Ks filter

# read in data
AO_lines = open('AOsensitivity.txt').readlines()
AO_sep = np.array([float(line.split(' ')[0]) for line in AO_lines[2:]]) # arcsec
AO_dmag = np.array([float(line.split(' ')[1]) for line in AO_lines[2:]]) # Ks filter

# convert separation to sma
AO_sma = np.array([arcsec2au(s,system_distance) for s in AO_sep])

# convert dmag to mass
AO_host_mag = mass2kmag(1.03)
AO_mag = AO_host_mag + AO_dmag
AO_mass = kmag2mass(AO_mag)



# ========= Speckle data ==========================================
# sep in arcsec, dmag @692nm = V band

# read in data
Speckle_lines = open('GeminiSpeckle692Sensitivity.txt').readlines()
Speckle_sep = np.array([float(line.split('  ')[0]) for line in Speckle_lines[12:]])
Speckle_dmag = np.array([float(line.split('  ')[1]) for line in Speckle_lines[12:]]) #692nm=Vband

# convert separation to sma
Speckle_sma = np.array([arcsec2au(s,186.9) for s in Speckle_sep])

# convert dmag to mass
Speckle_host_mag = mass2vmag(1.03)
Speckle_mag = Speckle_host_mag + Speckle_dmag
Speckle_mass = vmag2mass(Speckle_mag)



# ============ Gaia data =======================================
GAIA_sep = arcsec2au(1,system_distance)
GAIA_mass = 1


# =============== RV data ====================================
# trend in m/s/d converted to sma, mass

trend_msd = 0.006 # m/s/day
trend = trend_msd*356 # m/s/yr
max_trend = trend*2 # upper limit on trend
RV_sma = np.logspace(-2,5,100)
RV_mass = (max_trend/6.57)*(RV_sma/5)**2.*(Mjup/Msun)



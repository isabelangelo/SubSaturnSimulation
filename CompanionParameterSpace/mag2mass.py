import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# ===================================================
"""
Functions for converting magnitude to mass:

    kmag2mass: input mag (Kmag) ; dtype=float,array 
               output mass (M_sun) ; dtype=ndarray
    vmag2mass: input mag (Vmag) ; dtype=float,array 
               output mass (M_sun) ; dtype=ndarray
    
The functions are generated using a 1d linear spline interpolation of data from
Pecaut and Mamajek (2013)
"""

table = pd.read_csv('PecautMamajek_table.csv', delim_whitespace=True)
table = table.replace('...',np.nan)

M_K = np.array([float(i) for i in table['M_K']])
Mv = np.array([float(i) for i in table['Mv']])
M_sun = np.array([float(i) for i in table['Msun']])

kmag2mass = interp1d(M_K, M_sun)
vmag2mass = interp1d(Mv, M_sun)


# ===================================================
"""
Functions for converting maass to magnitude:

    mass2kmag: input mass (M_sun) ; dtype=float,array 
               output mag (Kmag) ; dtype=ndarray
    mass2vmag: input mass (M_sun) ; dtype=float,array 
               output mag (Vmag) ; dtype=ndarray
    
The functions are generated using a 1d linear spline interpolation of data from
Pecaut and Mamajek (2013)
"""

mass2kmag = interp1d(M_sun, M_K)
mass2vmag = interp1d(M_sun, Mv)

# ===================================================
from cgs import *

# ===================================================
def arcsec2au(sep, d):
    """
    Function to convert angular separation to physical distance
    :param sep: separation in arcsec
    :param d: Distance in parsecs
    :return: a in AU
    """
    
    sep_rad = sep*pi/(180*3600.)
    distance_cm = d*pc
    a_cm = distance_cm*sep_rad
    a = a_cm/AU
    return a
    
# ===================================================

def au2arcsec(a, d):
    """
    Function to convert angular separation to physical distance
    :param a: semi-major axis in AU
    :param d: Distance in parsecs
    :return: sep in arcsec
    """
    a_cm = a*AU
    distance_cm = d*pc
    sep_rad = a_cm/distance_cm
    sep_arcsec = sep_rad*(180*3600)/pi
    return sep_arcsec
    
# ===================================================
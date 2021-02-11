# Created a custom config, because default did not converge on high eccentricity

import radvel
import numpy as np
import pandas as pd
#from rvsearch import utils

starname = 'CK00367'
nplanets = 2
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', 2: 'c',}

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc secosw sesinw k', 
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=31.5786453, vary=False)
anybasis_params['tc1'] = radvel.Parameter(value=2457820.52, vary=False)
anybasis_params['secosw1'] = radvel.Parameter(value=0.5810851924490481, vary=False)
anybasis_params['sesinw1'] = radvel.Parameter(value=0.6901738904917025, vary=False)
anybasis_params['k1'] = radvel.Parameter(value=16.7, vary=False)

# Define prior centers for second planet
anybasis_params['per2'] = radvel.Parameter(value=1632.0)
anybasis_params['tc2'] = radvel.Parameter(value=2458922)
anybasis_params['secosw2'] = radvel.Parameter(value=0.00)
anybasis_params['sesinw2'] = radvel.Parameter(value=0.00)
anybasis_params['k2'] = radvel.Parameter(value=1)

time_base = 2457806.751140
anybasis_params['dvdt'] = radvel.Parameter(value=0.0,vary=False)
anybasis_params['curv'] = radvel.Parameter(value=0.0,vary=False)

# read in data
data = pd.read_csv('CK00367_default_data.csv',
                   dtype={'time': np.float64, 'mnvel': np.float64, 'err': np.float64, 'tel': str})
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.05)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel

# store telescope information
instnames = ['j']
ntels = len(instnames)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=True, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)

# initialize model
params = anybasis_params.basis.to_any_basis(anybasis_params,fitting_basis)
mod = radvel.RVModel(params, time_base=time_base)

# add priors for planets
priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          radvel.prior.HardBounds('jit_j', 0.0, 10.0),
          radvel.prior.Gaussian('per1', 1600, 400)
          ]

stellar = dict(mstar=1.03, mstar_err=0.04)


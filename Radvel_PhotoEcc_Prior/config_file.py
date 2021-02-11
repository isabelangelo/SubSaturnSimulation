# Created a custom config, because default did not converge on high eccentricity

import radvel
import numpy as np
import pandas as pd
#from rvsearch import utils

starname = 'CK00367'
nplanets = 1
fitting_basis = 'per tc secosw sesinw k'
bjd0 = 0.
planet_letters = {1: 'b', }

# Define prior centers (initial guesses) in a basis of your choice (need not be in the fitting basis)
anybasis_params = radvel.Parameters(nplanets, basis='per tc e w k', 
                                    planet_letters=planet_letters) # initialize Parameters object
anybasis_params['per1'] = radvel.Parameter(value=31.578659)
anybasis_params['tc1'] = radvel.Parameter(value=2457820.706751)
anybasis_params['e1'] = radvel.Parameter(value=0.836000)
anybasis_params['w1'] = radvel.Parameter(value=0.940732)
anybasis_params['k1'] = radvel.Parameter(value=17.600000)

time_base = 2457806.751140
anybasis_params['dvdt'] = radvel.Parameter(value=0.0)
anybasis_params['curv'] = radvel.Parameter(value=0.0)
data = pd.read_csv('CK00367_data.csv',
                   dtype={'time': np.float64, 'mnvel': np.float64, 'err': np.float64, 'tel': str})
bin_t, bin_vel, bin_err, bin_tel = radvel.utils.bintels(data['time'].values, data['mnvel'].values, data['errvel'].values, data['tel'].values, binsize=0.05)
data = pd.DataFrame([], columns=['time', 'mnvel', 'errvel', 'tel'])
data['time'] = bin_t
data['mnvel'] = bin_vel
data['errvel'] = bin_err
data['tel'] = bin_tel

instnames = ['j']
ntels = len(instnames)
anybasis_params['gamma_j'] = radvel.Parameter(value=0.0, vary=True, linear=True)
anybasis_params['jit_j'] = radvel.Parameter(value=1.0)

params = anybasis_params.basis.to_any_basis(anybasis_params,fitting_basis)
mod = radvel.RVModel(params, time_base=time_base)

mod.params['per1'].vary = True
mod.params['tc1'].vary = True
mod.params['secosw1'].vary = True
mod.params['sesinw1'].vary = True
mod.params['dvdt'].vary = True
mod.params['curv'].vary = False
mod.params['jit_j'].vary = True

# add rho posterior from transits as prior
rho_post = pd.read_csv('./' + 'CK00367' + planet_letters[1] + '-rho_ratio_interp-ga_plx_teff_feh.csv')

priors = [
          radvel.prior.EccentricityPrior(nplanets),
          radvel.prior.PositiveKPrior(nplanets),
          radvel.prior.HardBounds('jit_j', 0.0, 10.0),
          radvel.prior.Gaussian('per1', 31.5786452580379, 2e-07), # where does this come from?
          radvel.prior.Gaussian('tc1', 2457820.706751, 0.5),
          radvel.prior.NumericalPrior_PhotoEcc(param_list=["secosw1", "sesinw1"], values=rho_post)
          ]

stellar = dict(mstar=1.03, mstar_err=0.04)



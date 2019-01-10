import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot()

import numpy as np
config.numpy()
from collections import OrderedDict as odict
from copy import deepcopy
import system

def get_sims_structure():
  # build a dictionary of model structural variants, starting from the most complicated
  specs = system.get_specs()
  sims = odict()
  for nu in [specs['params']['nu'],specs['params']['mu']]:
    for G in [3,1]:
      for Z in [1,0]:
        # get the default model
        model = system.get_model()
        # define the growth rate
        model.params['nu'].update(nu)
        # define the number of of groups
        if G == 1:
          if Z != 0: continue # Z is irrelevant for G = 1, so only need 1 variant
          model.collapse(['ii'])
        # define turnover via specified durations
        if Z == 0:
          model.params['dur'].update(np.nan)
          model.params['zeta'].update(np.nan)
        # add simulation to dictionary
        sims.update({
          'nu={}_G={}_Z={}'.format(nu,G,Z):
          system.get_simulation(model)
        })
  return sims

def get_sims_zeta():
  # build a dictionary of model turnover variants, starting with none

  def get_zeta_specs(nu,mu,pe,dur,zeta):
    # build set of sufficient constraints given zeta with some entries zero:
    #   zeta will be either zero or nan (calculated)
    #   dur will be either input value of nan (calculated)
    #   pe will be nan (calculated) always
    # NOTE did not check validity of this constraint building for G > 3
    zeta = deepcopy(zeta)
    for i,zi in enumerate(zeta):
      zc = sum(np.isnan(zi))
      if zc == 0:
        # calculate dur if all zeta for this i specified
        dur[i] = np.nan
      if zc == 2:
        # specify zeta for this i if all to be calculated (even split)
        zsum = (1/dur[i] - mu)
        zeta[i] = [zsum/2 if (i != j) else (0) for j in range(3)]
      # calculate pe always
      pe[i] = np.nan
    return {'zeta':zeta,'dur':dur,'pe':pe}

  def zeta_str(zeta):
    # make a string to represent the variant like [[.ZZ][0.*][00.]]
    return str([[
        '.' if i==j else
        '0' if zij == 0 else
        'x' if np.isnan(zij) else
        'Z'
        for j,zij in enumerate(zi)]
      for i,zi in enumerate(zeta)]
    ).replace('\'','').replace(', ','')
  
  specs = system.get_specs()
  sims = odict()
  zeta = np.zeros((3,3))
  # adding turnover flows one at a time
  for i,ij in enumerate([None,(0,2),(0,1),(1,2),(2,0),(2,1),(1,0)]):
    model = system.get_model()
    if ij:
      zeta[ij[0],ij[1]] = np.nan
    # define the sufficient constraints for this system
    for ki in ['M','W']:
      zspecs = get_zeta_specs(nu = model.params['nu'],
                              mu = model.params['mu'],
                              pe = model.params['pe'].islice(ki=ki),
                              dur = model.params['dur'].islice(ki=ki),
                              zeta = zeta)
      for param in ['zeta','dur','pe']:
        model.params[param].update(zspecs[param],ki=ki)
    # add simulation to dictionary
    sims.update({
      '({})-{}'.format(i,zeta_str(zspecs['zeta'])):
      system.get_simulation(model)
    })
  return sims

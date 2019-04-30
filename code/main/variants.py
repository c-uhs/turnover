import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
config.epimodel()
config.plot()

import numpy as np
config.numpy()
from collections import OrderedDict as odict
from copy import deepcopy
import re
import system
import elements

mu = system.get_specs()['params']['mu']

def iter_params():
  specs = system.get_specs()
  for nu in [specs['params']['nu'],specs['params']['mu']]:
    for G in [3.,1.]:
      for Z in [1.,0.]:
        if not ((G == 1) and (Z == 1)):
          yield nu,G,Z

def make_name(nu,G,Z):
  return 'nu={}_G={}_Z={}'.format(float(nu),int(G),int(Z))

def parse_name(name):
  nu,G,Z = re.findall('nu\=(.*)\_G\=(.*)\_Z\=(.*)',name)[0]
  return float(nu),int(G),int(Z)

def make_title(nu,G,Z,i=None):
  return '{}$\\nu {} \\mu, G = {}, \\zeta {} 0$'.format(
    'V{}: '.format(i+1) if (i is not None) else '',
    '>' if (nu > mu) else '=',
    int(G),
    '>' if bool(Z) else '=',
  )

def get_sims(t=None):
  # build a dictionary of model structural variants, starting from the most complicated
  specs = system.get_specs()
  sims = odict()
  for nu,G,Z in iter_params():
    # get the default model
    model = system.get_model()
    # define the growth rate
    model.params['nu'].update(nu)
    # define the number of of groups
    if G == 1:
      model.collapse(['ii'])
    # define turnover via specified durations
    if Z == 0:
      model.params['dur'].update(np.nan)
      model.params['zeta'].update(np.nan)
    # add simulation to dictionary
    sims.update({
      make_name(nu,G,Z):
      system.get_simulation(model,t=t)
    })
  return sims

def get_selectors():
  select = odict()
  for nu,G,Z in iter_params():
    name = make_name(nu,G,Z)
    if nu == list(iter_params())[0][0]:
      color = elements.Color([0.8,  0,  0]).lighten(0.6*(1-Z))
    else:
      color = elements.Color([  0,  0,0.8]).lighten(0.6*(1-Z))
    linestyle = '--' if (G == 1) else '-'
    select.update({
      name:
      elements.Selector(
        name = name,
        select = {},
        title = make_title(nu,G,Z),
        color = color,
        linestyle = linestyle,
      )
    })
  return select

def iter_all(t=None):
  select = get_selectors()
  sims = get_sims(t=t)
  for name,sim in sims.items():
    yield name,sim,select[name]

import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
import numpy as np
config.epimodel()
config.plot(tex=True)
import matplotlib.pyplot as plt
from collections import OrderedDict as odict
import modelutils
import elements
import system
import utils

def iter_dh(N,dlims=None):
  dlims = dlims if dlims is not None else [33,3]
  for dh in np.logspace(np.log10(dlims[0]),np.log10(dlims[1]),N):
    yield np.around(dh,3)

def iter_tau(N,tlims=None):
  tlims = tlims if tlims is not None else [0.0,0.2]
  for tau in np.linspace(tlims[0],tlims[1],N):
    yield np.around(tau,3)

def iter_2d(N,dlims=None,tlims=None):
  for it,tau in enumerate(iter_tau(N,tlims)):
    for id,dh in enumerate(iter_dh(N,dlims)):
      yield id,it,dh,tau

def get_selectors(sim):
  names = ['all','high','med','low']
  for name in names:
    yield sim.model.select[name]

def get_sim(dh,tau):
  model = system.get_model()
  dhmax = 1.0/model.params['mu']
  model.params['dur'].update([dh,min(dh*5,dhmax),min(dh*30,dhmax)])
  z1 = (1/model.params['dur'].iselect(ii='H',ki='M') - model.params['mu'])/2
  z2 = (1/model.params['dur'].iselect(ii='M',ki='M') - model.params['mu'])/2
  z3 = (1/model.params['dur'].iselect(ii='L',ki='M') - model.params['mu'])/2
  model.params['zeta'].update(z1,ii='H',ip='M')
  model.params['zeta'].update(z1,ii='H',ip='L')
  model.params['zeta'].update(z2,ii='M',ip='L')
  model.params['zeta'].update(z2,ii='M',ip='H')
  model.params['zeta'].update(z3,ii='L',ip='M')
  model.params['zeta'].update(z3,ii='L',ip='H')
  # model.params['zeta'].update(np.nan)
  model.params['pe'].update(np.nan)
  model.params['tau'].update(tau)
  dt = min(dh*0.8,0.2)
  t = system.get_t(dt=dt,tmin=0,tmax=200)
  sim = system.get_simulation(model,t=t)
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = ['prevalence','incidence'],
  ))
  return sim

def fname_data(output,select,N):
  return os.path.join(config.path['specs'],'2d','N={}'.format(N),
    '{}-{}.csv'.format(output,select))

def run_sim(sim,dh,tau,id,it,N):
  def save(output,select):
    fname = fname_data(output.name,select.name,N)
    if not os.path.exists(fname):
      data = [[None for x in range(N)] for y in range(N)]
    else:
      data = utils.loadcsv(fname,asdict=False)
    data[it][id] = modelutils.taccum(output,**select)[0]
    utils.savecsv(fname,data,append=False)
  sim.solve()
  tmax = sim.t[-1]
  for output in sim.outputs.values():
    for select in get_selectors(sim):
      select.update(t=tmax)
      save(output,select)

if __name__ == '__main__':
  N = 15
  idx = int(sys.argv[1])
  for id,it,dh,tau in iter_2d(N):
    if id == idx:
      print((id,it,dh,tau),flush=True)
      sim = get_sim(dh,tau)
      run_sim(sim,dh,tau,id,it,N)

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

N = 15
OUTPUTS = ['prevalence','incidence']
SELECTORS = ['all','high','med','low']

# iteration functions

def iter_dh(N,dlims=None):
  dlims = dlims if dlims is not None else [33,3]
  for dh in np.logspace(np.log10(dlims[0]),np.log10(dlims[1]),N):
    yield np.around(dh,3)

def iter_tau(N,tlims=None):
  tlims = tlims if tlims is not None else [0.0,0.2]
  for tau in np.linspace(tlims[0],tlims[1],N):
    yield np.around(tau,3)

def iter_di(N,tlims=None):
  for tau in iter_tau(N,tlims):
    yield np.around(1/tau,1)

def iter_both(N,dlims=None,tlims=None):
  for it,tau in enumerate(iter_tau(N,tlims)):
    for id,dh in enumerate(iter_dh(N,dlims)):
      yield id,it,dh,tau

def iter_selectors(sim):
  for name in SELECTORS:
    yield sim.model.select[name]

# filenames for saving (data, figures)

def fname_data(output,select):
  return os.path.join(config.path['data'],'surface',
    '{}-{}.csv'.format(output,select))

def fname_fig(output,select,norm):
  return os.path.join(config.path['figs'],'plots','surface',
    '2d-{}{}-{}.eps'.format(output,'-norm' if norm else '',select))

# define the simulation parameters

def get_sim(dh,tau):
  model = system.get_model()
  dhmax = 1.0/model.params['mu']
  model.params['dur'].update([dh,min(dh*5,dhmax),min(dh*30,dhmax)])
  z1 = (1/model.params['dur'].iselect(ii=['H'],ki=['M']) - model.params['mu'])/2
  z2 = (1/model.params['dur'].iselect(ii=['M'],ki=['M']) - model.params['mu'])/2
  z3 = (1/model.params['dur'].iselect(ii=['L'],ki=['M']) - model.params['mu'])/2
  model.params['zeta'].update(z1,ii=['H'],ip=['M'])
  model.params['zeta'].update(z1,ii=['H'],ip=['L'])
  model.params['zeta'].update(z2,ii=['M'],ip=['L'])
  model.params['zeta'].update(z2,ii=['M'],ip=['H'])
  model.params['zeta'].update(z3,ii=['L'],ip=['M'])
  model.params['zeta'].update(z3,ii=['L'],ip=['H'])
  model.params['pe'].update(np.nan)
  model.params['tau'].update(tau)
  dt = min(dh*0.8,0.2)
  t = system.get_t(dt=dt,tmin=0,tmax=200)
  sim = system.get_simulation(model,t=t)
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = OUTPUTS,
  ))
  return sim

# run sims and save the 'steady-state' outputs

def run_sims(idx=[]):
  idx = [int(i) for i in idx]
  for id,it,dh,tau in iter_both(N):
    if id in idx:
      print('Dh = {:6.3f} | tau = {:6.3f}'.format(dh,tau),flush=True)
      sim = get_sim(dh,tau)
      run_sim(sim,id,it,N)

def run_sim(sim,id,it,N):
  def save(output,select):
    fname = fname_data(output.name,select.name)
    if not os.path.exists(fname):
      data = [[None for x in range(N)] for y in range(N)]
    else:
      data = utils.loadcsv(fname,asdict=False)
    data[it][id] = modelutils.taccum(output,**select)[0]
    utils.savecsv(fname,data,append=False)
  sim.solve()
  tmax = sim.t[-1]
  for output in sim.outputs.values():
    for select in iter_selectors(sim):
      select.update(t=[tmax])
      save(output,select)

# make surface plots of the result

def make_plot(output,select,N,norm=False):
  def clean_data(data,norm=False):
    return [[
      np.float(d) / (1 if not norm else np.float(da[0]))
      for d in da]
      for da in data]
  def ticks(values,incr,dec):
    return range(len(values))[::incr],[np.around(v,dec) for v in values][::incr]
  data = clean_data(utils.loadcsv(fname_data(output,select),asdict=False),norm)
  plt.figure(figsize=(6.5,5))
  plt.imshow(data,cmap=plt.get_cmap('inferno'),interpolation='none')
  plt.colorbar().ax.tick_params(labelsize=16)
  plt.yticks(*ticks(list(iter_di(N)),2,2),fontsize=16)
  plt.xticks(*ticks(list(iter_dh(N)),2,1),fontsize=16)
  plt.ylabel('Duration of Infectiousness $\\delta_{\\mathcal{I}}$',fontsize=20)
  plt.xlabel('Duration in High-Risk Group $\\delta_H$',fontsize=20)
  plt.grid(None)
  plt.tight_layout(pad=0.5)
  plt.text(N-1,N-1,'$R_0 < 1$',fontsize=20,color='w',va='bottom',ha='right')
  plt.savefig(fname_fig(output,select,norm))
  # plt.show()
  plt.close()

def make_plots():
  for output in OUTPUTS:
    for select in SELECTORS:
      make_plot(output,select,N,norm=False)
      make_plot(output,select,N,norm=True)

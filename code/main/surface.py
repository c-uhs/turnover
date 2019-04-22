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

N = 31
OUTPUTS = ['prevalence','incidence']
SELECTORS = ['all','high','med','low']

# iteration functions

def iter_dh(N,dhlims=None):
  dhlims = dhlims if dhlims is not None else [3,33]
  for dur in np.logspace(*np.log10(dhlims),N):
    yield np.around(dur,3)

def iter_di(N,dilims=None):
  dilims = dilims if dilims is not None else [1,20]
  for dur in np.linspace(*dilims,N)[::-1]:
    yield np.around(dur,3)

def iter_both(N,dhlims=None,dilims=None):
  for idi,di in enumerate(iter_di(N,dilims)):
    for idh,dh in enumerate(iter_dh(N,dhlims)):
      yield idh,idi,dh,di

def iter_selectors(sim):
  for name in SELECTORS:
    yield sim.model.select[name]

# filenames for saving (data, figures)

def fname_data(output,select,dh,di):
  return os.path.join(config.path['data'],'surface',output,select,
    'dh={}-di={}.txt'.format(dh,di))

def fname_fig(output,select,norm):
  return os.path.join(config.path['figs'],'plots','surface',
    'surface-{}{}-{}.eps'.format(output,'-norm' if norm else '',select))

# define the simulation parameters

def get_sim(dh,di,tmax=200):
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
  model.params['tau'].update(1/di)
  dt = min(dh*0.8,0.5)
  t = system.get_t(dt=dt,tmin=0,tmax=tmax)
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
  for idh,idi,dh,di in iter_both(N):
    if idi in idx:
      print('Dh = {:6.3f} | Di = {:6.3f}'.format(dh,di),flush=True)
      sim = get_sim(dh,di)
      run_sim(sim,dh,di,N)

def run_sim(sim,dh,di,N):
  sim.solve()
  tmax = sim.t[-1]
  for output in sim.outputs.values():
    for select in iter_selectors(sim):
      select.update(t=[tmax])
      fname = fname_data(output.name,select.name,dh,di)
      value = modelutils.taccum(output,**select)[0]
      utils.savecsv(fname,[[value]],append=False)

# make surface plots of the result

def load_data(output,select,N,norm=False):
  dhn = list(iter_dh(N))[-1]
  def load_one(dh,di):
    x = np.float(utils.loadcsv(fname_data(output,select,dh,di),asdict=False)[0][0])
    n = np.float(utils.loadcsv(fname_data(output,select,dhn,di),asdict=False)[0][0]) if norm else 1
    return (x / n) if x > 1e-5 else 0
  return [[ load_one(dh,di) for dh in iter_dh(N)] for di in iter_di(N)]

def make_surface(data,N,labels=None,fs=16):
  def ticks(values,incr,dec):
    return range(len(values))[::incr],[np.around(v,dec) for v in values][::incr]
  cmap = plt.get_cmap('inferno')
  plt.imshow(data,cmap=cmap,interpolation='none')
  plt.colorbar().ax.tick_params(labelsize=fs)
  plt.yticks(*ticks(list(iter_di(N)),5,2),fontsize=fs)
  plt.xticks(*ticks(list(iter_dh(N)),5,1),fontsize=fs)
  if labels in [None,'long']:
    plt.ylabel('Duration of Infectiousness $\\delta_{\\mathcal{I}}$',fontsize=fs)
    plt.xlabel('Duration in High-Risk Group $\\delta_H$',fontsize=fs)
  elif labels in ['short']:
    plt.ylabel('$\\delta_{\\mathcal{I}}$')
    plt.xlabel('$\\delta_H$')
  plt.grid(None)

def make_plot(data,N,save=None):
  plt.figure(figsize=(6.5,5))
  make_surface(data,N)
  plt.tight_layout(pad=0.5)
  plt.text(0,N-1,'$R_0 < 1$',fontsize=20,color='w',va='bottom',ha='left')
  if save:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def make_plots():
  for output in OUTPUTS:
    for select in SELECTORS:
      for norm in [False,True]:
        make_plot(load_data(output,select,N,norm),N,fname_fig(output,select,norm))

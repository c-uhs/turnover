import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
import numpy as np
config.epimodel()
config.plot(tex=True)
import matplotlib.pyplot as plt
from scipy import interpolate
from collections import OrderedDict as odict
import modelutils
import elements
import system
import utils

N = 31
OUTPUTS = ['prevalence','incidence']
SELECTORS = ['all','high','med','low']
TAU1D = 0.1
CA = 3 # HACK
SAVE = False

# iteration functions

def iter_phi(N,plims=None):
  plims = plims if plims is not None else [1/3,0.03]
  for phi in 1/np.logspace(*np.log10(1/np.array(plims)),N)[::-1]:
    yield np.around(phi,3)

def iter_tau(N,tlims=None):
  tlims = tlims if tlims is not None else [0.05,1]
  for tau in 1/np.linspace(*1/np.array(tlims),N)[::-1]:
    yield np.around(tau,3)

def iter_both(N,plims=None,tlims=None):
  for it,tau in enumerate(iter_tau(N,tlims)):
    for ip,phi in enumerate(iter_phi(N,plims)):
      yield ip,it,phi,tau

def iter_selectors(sim):
  for name in SELECTORS:
    yield sim.model.select[name]

# filenames for saving (data, figures)

def fname_fun(*args,**kwargs):
  args = list(args)
  keys = sorted(list(kwargs.keys()))
  for key in keys:
    value = kwargs[key]
    if isinstance(value,bool):
      if value:
        args.append(key)
    else:
      args.append('{}={}'.format(key,value))
  return '-'.join(args)

def fname_element(output,select,*args,**kwargs): # output,select,phi=None,tau=None
  return os.path.join(config.path['data'],'surface','raw',output,select,
    fname_fun(*args,**kwargs)+'.txt')

def fname_array(*args,**kwargs): # output,select,norm=False,tau=None
  return os.path.join(config.path['data'],'surface',
    fname_fun(*args,**kwargs)+'.csv')

def fname_fig(*args,**kwargs):
  return os.path.join(config.path['figs'],'plots','sensitivity',
    fname_fun(*args,**kwargs)+'.pdf')

# utils

def ticks(values,incr,dec):
  return list(range(len(values)))[::incr], [np.around(v,dec) for v in values][::incr]

# define the simulation parameters

def get_sim(phi,tau,tmax=200):
  model = system.get_model()
  dmax = 1.0/model.params['mu']
  dh = 1/phi
  model.params['dur'].update(dh,               ii=['H'])
  model.params['dur'].update(dh+0.30*(dmax-dh),ii=['M'])
  model.params['dur'].update(np.nan,           ii=['L'])
  z1 = (1/model.params['dur'].iselect(ii=['H'],ki=['M']) - model.params['mu'])/2
  model.params['phi'].update(z1,ii=['H'],ip=['M'])
  model.params['phi'].update(z1,ii=['H'],ip=['L'])
  model.params['tau'].update(tau)
  dt = min(dh*0.8,0.5)
  t = system.get_t(dt=dt,tmin=0,tmax=tmax)
  sim = system.get_simulation(model,t=t)
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = OUTPUTS+['C','X'],
  ))
  return sim

# run sims and save the equilibrium outputs

def run_sims(idx=[]):
  idx = [int(i) for i in idx]
  for ip,it,phi,tau in iter_both(N):
    if it in idx:
      print('phi = {:6.3f} | tau = {:6.3f}'.format(phi,tau),flush=True)
      sim = get_sim(phi,tau)
      run_sim(sim,phi,tau)
  # HACK: run tau = 0.1 results
  # tau = 0.1
  # for phi in iter_phi(N):
  #   sim = get_sim(phi,tau)
  #   run_sim(sim,phi,tau)

def run_sim(sim,phi,tau):
  sim.solve()
  tmax = sim.t[-1]
  for output in OUTPUTS:
    for select in iter_selectors(sim):
      save_element(sim.outputs[output],select,tmax,phi,tau)
  for output in ['C','X']:
    save_element(sim.outputs[output],sim.model.select['I'],tmax,phi,tau)

def save_element(output,select,tmax,phi,tau):
  fname = fname_element(output.name,select.name,phi=phi,tau=tau)
  value = modelutils.taccum(output,**select).islice(t=tmax)
  utils.savetxt(fname,float(value))

# make surface plots of the result

def load_element(output,select,phi,tau,norm):
  phin = list(iter_phi(N))[0]
  x = np.float(utils.loadtxt(fname_element(output,select,phi=phi,tau=tau)))
  n = np.float(utils.loadtxt(fname_element(output,select,phi=phin,tau=tau))) if norm else 1
  return (x / n) if x > 1e-5 else 0

def load_data(output,select,norm=False,phi=None,tau=None):
  phi = phi if phi is not None else list(iter_phi(N))
  tau = tau if tau is not None else list(iter_tau(N))
  dname = fname_array(output,select,norm=norm)
  if True: # not os.path.exists(dname):
    array = [[ load_element(output,select,p,t,norm) for p in phi] for t in tau]
    utils.savecsv(dname,array,append=False)
  else:
    array = utils.loadcsv(dname,asdict=False)
  array = np.array(array,dtype=np.float)
  return array

def make_1d_pretty(output):
  plt.xticks(*ticks(list(iter_phi(N)),5,2))
  plt.xlabel('High-Risk Turnover $\\delta_H^{-1}$')
  ylabel = {
    'prevalence': 'Prevalence',
    'incidence':  'Incidence (per 1000 person-years)',
    'X':          'Proportion of population who are infectious',
    'C':          'Average contact rate of infectious individuals',
    'XC':         'Infectious partnerships proportion',
  }
  plt.ylabel(ylabel[output])

def draw_regions(xp,yp,labels):
  yl = plt.ylim()
  xl = plt.xlim()
  x  = [xp[0]] + xp + [xp[-1], xl[1], xl[1]]
  y  = [yl[0]] + yp + [yl[ 1], yl[1], yl[0]]
  plt.fill(x, y, color = 'r', alpha = 0.1)
  lp = lambda lim,p: lim[0]+p*(lim[1]-lim[0])
  tprops = dict(va='top',bbox=dict(boxstyle='round',fc='w',ec='0.8',alpha=0.9))
  plt.text(lp(xl,0.05),lp(yl,0.95),labels[0],ha='left',**tprops)
  plt.text(lp(xl,0.95),lp(yl,0.95),labels[1],ha='right',**tprops)
  plt.xlim(xl)
  plt.ylim(yl)

def make_1d(data,output,regions=False):
  def get_peaks(data,interp):
    x  = np.arange(0,N,1)
    xi = np.arange(0,N,interp)
    fi = lambda y: interpolate.interp1d(x,y,kind='cubic',bounds_error=False)(xi)
    xp = [np.nanargmax(fi(y))*interp for y in data]
    yp = [np.nanmax(fi(y)) for y in data]
    return xp,yp
  plt.gca().set_prop_cycle('color',plt.cm.Blues_r(np.linspace(0,1,data.shape[0])))
  plt.plot(range(0,N),data.transpose())
  if regions:
    xp,yp = get_peaks(data,0.01)
    if np.any(xp):
      draw_regions(xp,yp,['A','B'])
  make_1d_pretty(output)

def make_surface(data,labels=None):
  cmap = plt.get_cmap('inferno')
  plt.imshow(data,cmap=cmap,interpolation='none')
  plt.colorbar()
  plt.yticks(*ticks(list(iter_tau(N)),5,2))
  plt.xticks(*ticks(list(iter_phi(N)),5,2))
  if labels in [None,'long']:
    plt.ylabel('Rate of Treatment $\\tau$')
    plt.xlabel('High-Risk Turnover $\\delta_H^{-1}$')
  elif labels in ['short']:
    plt.ylabel('$\\tau$')
    plt.xlabel('$\\phi$')
  plt.grid(None)

def gen_2d_plot(data,save=None):
  plt.figure(figsize=(4,3))
  make_surface(data)
  plt.tight_layout(pad=0.2)
  plt.text(N-1,1,'$R_0 < 1$',fontsize=14,color='w',va='top',ha='right')
  if save and SAVE:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def gen_1d_plot(data,output,save=None,regions=False):
  plt.figure(figsize=(4,3))
  make_1d(data,output,regions=regions)
  plt.gca().set_position([0.16,0.14,0.82,0.84])
  if save and SAVE:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def gen_plots():
  for output in OUTPUTS:
    for select in ['low']:
      gen_1d_plot(load_data(output,select,tau=[TAU1D]),output,fname_fig('1d',output,select,tau=TAU1D),regions=True)
      for norm in [False,True]:
        gen_2d_plot(load_2d_data(output,select,norm),fname_fig('surface',output,select,norm=norm))
  data = {
    '1d': {output: load_data(output,'I',tau=[TAU1D]) for output in ['X','C'] },
    '2d': {output: load_data(output,'I')             for output in ['X','C'] },
  }
  data['1d']['XC'] = data['1d']['X'] * data['1d']['C'] / CA
  data['2d']['XC'] = data['2d']['X'] * data['2d']['C'] / CA
  for output in ['X','C','XC']:
    gen_1d_plot(data['1d'][output],output,fname_fig('1d',output,'I',tau=TAU1D),regions=True)
    gen_1d_plot(data['2d'][output],output,fname_fig('2d',output,'I'),regions=True)
  gen_1d_plot(load_data('incidence','all'),'incidence',fname_fig('2d','incidence','all'),regions=True)

import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
import numpy as np
config.epimodel()
config.plot(tex=True)
import matplotlib.pyplot as plt
from scipy import interpolate
from collections import OrderedDict as odict
from textwrap import wrap
import modelutils
import elements
import system
import utils

OUTPUTS = ['prevalence','incidence']
SELECTORS = ['all','high','med','low']
TAU1D = 0.1
CA = 3 # HACK

# iteration functions

def iter_phi(plims=None):
  plims = plims if plims is not None else [1/3,0.03]
  for phi in 1/np.logspace(*np.log10(1/np.array(plims)),config.N)[::-1]:
    yield np.around(phi,3)

def iter_tau(tlims=None):
  tlims = tlims if tlims is not None else [0.05,1]
  for tau in 1/np.linspace(*1/np.array(tlims),config.N)[::-1]:
    yield np.around(tau,3)

def iter_both(plims=None,tlims=None):
  for it,tau in enumerate(iter_tau(tlims)):
    for ip,phi in enumerate(iter_phi(plims)):
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

def fname_element(output,select,*args,**kwargs):
  return os.path.join(config.path['data'],'sensitivity','raw',output,select,
    fname_fun(*args,**kwargs)+'.txt')

def fname_array(*args,**kwargs):
  return os.path.join(config.path['data'],'sensitivity',
    fname_fun(*args,**kwargs)+'.csv')

def fname_fig(*args,**kwargs):
  return os.path.join(config.path['figs'],'sensitivity',
    fname_fun(*args,**kwargs)+'.pdf')

# utils

def ticks(values,incr,dec):
  vfun = lambda v: float(v) if dec else int(v)
  return list(range(len(values)))[::incr], [np.around(vfun(v),dec) for v in values][::incr]

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
    names = OUTPUTS+['C','X','tip'],
  ))
  return sim

# run sims and save the equilibrium outputs

def run_sims(idx=-1):
  if int(idx) == -1:
    for phi in iter_phi():
      run_sim(phi=phi,tau=TAU1D)
  else:
    for ip,it,phi,tau in iter_both():
      if it == int(idx):
        run_sim(phi=phi,tau=tau)

def run_sim(phi,tau):
  print('phi = {:6.3f} | tau = {:6.3f}'.format(phi,tau),flush=True)
  sim = get_sim(phi,tau)
  sim.solve()
  tmax = sim.t[-1]
  for output in OUTPUTS:
    for select in iter_selectors(sim):
      save_element(sim.outputs[output],select,tmax,phi,tau)
  for output in ['C','X']:
    save_element(sim.outputs[output],sim.model.select['I'],tmax,phi,tau)
  for select in iter_selectors(sim):
    save_element(sim.outputs['tip'],select,tmax,phi,tau)

def save_element(output,select,tmax,phi,tau):
  fname = fname_element(output.name,select.name,phi=phi,tau=tau)
  value = modelutils.taccum(output,**select).islice(t=tmax)
  utils.savetxt(fname,float(value))

# make surface plots of the result

def load_element(output,select,phi,tau,norm):
  phin = list(iter_phi())[0]
  x = np.float(utils.loadtxt(fname_element(output,select,phi=phi,tau=tau)))
  n = np.float(utils.loadtxt(fname_element(output,select,phi=phin,tau=tau))) if norm else 1
  return (x / n) if x > 1e-5 else 0

def load_data(output,select,norm=False,phi=None,tau=None):
  phi = phi if phi is not None else list(iter_phi())
  tau = tau if tau is not None else list(iter_tau())
  dname = fname_array(output,select,norm=norm)
  if True: # not os.path.exists(dname): # TODO: whats going on here?
    array = [[ load_element(output,select,p,t,norm) for p in phi] for t in tau]
    utils.savecsv(dname,array,append=False)
  else:
    array = utils.loadcsv(dname,asdict=False)
  array = np.array(array,dtype=np.float)
  return array

def make_1d_pretty(output,select=None):
  ytitle = {
    'prevalence': 'Prevalence',
    'incidence':  'Incidence (per 1000 PY)',
    'X':          'Proportion of population who are infectious',
    'C':          'Average contact rate of infectious individuals',
    'XC':         'Infectious partnerships proportion',
    'tip':        'Proportion of new infections from turnover',
  }
  yspec = {
    'high': ' among High-Risk',
    'med':  ' among Medium-Risk',
    'low':  ' among Low-Risk',
    'all':  ' Overall',
    None: '',
  }
  ylabel = ytitle[output]+yspec[select]
  if config.context == 'paper':
    plt.xticks(*ticks(list(iter_phi()),5,2))
    plt.xlabel('High-Risk Turnover $\\delta_H^{-1}$')
    plt.ylabel(ylabel)
  if config.context == 'isstdr':
    plt.xticks(*ticks([1/phi for phi in iter_phi()],5,0))
    plt.xlabel('Duration in High-Risk Group',
      fontsize='x-large',
      labelpad=5,
    )
    plt.ylabel('\n'.join(wrap(ylabel,13,break_on_hyphens=False)),
      fontsize='x-large',
      labelpad=35,
      rotation=0,
      va='center',
    )

def draw_regions(xp,yp,labels):
  yl = plt.ylim()
  xl = plt.xlim()
  x  = [xp[0]] + xp + [xp[-1], xl[1], xl[1]]
  y  = [yl[0]] + yp + [yl[ 1], yl[1], yl[0]]
  plt.fill(x, y, color = [1.0,0.0,0.0], alpha = 0.1)
  lp = lambda lim,p: lim[0]+p*(lim[1]-lim[0])
  tprops = dict(va='top',bbox=dict(boxstyle='round',fc='w',ec='0.8',alpha=0.9))
  plt.text(lp(xl,0.05),lp(yl,0.95),labels[0],ha='left',**tprops)
  plt.text(lp(xl,0.95),lp(yl,0.95),labels[1],ha='right',**tprops)
  plt.xlim(xl)
  plt.ylim(yl)

def make_1d(data,output,select,regions=False):
  def get_peaks(data,interp):
    x  = np.arange(0,config.N,1)
    xi = np.arange(0,config.N,interp)
    fi = lambda y: interpolate.interp1d(x,y,kind='cubic',bounds_error=False)(xi)
    xp = [np.nanargmax(fi(y))*interp for y in data]
    yp = [np.nanmax(fi(y)) for y in data]
    return xp,yp
  plt.gca().set_prop_cycle('color',plt.cm.Blues_r(np.linspace(0.1,0.9,data.shape[0])))
  plt.plot(range(0,config.N),data.transpose())
  if regions:
    xp,yp = get_peaks(data,0.01)
    if np.any(xp):
      draw_regions(xp,yp,['A','B'])
  make_1d_pretty(output,select)

def make_surface(data,labels=None):
  cmap = plt.get_cmap('inferno')
  plt.imshow(data,cmap=cmap,interpolation='none')
  plt.colorbar()
  plt.yticks(*ticks(list(iter_tau()),5,2))
  plt.xticks(*ticks(list(iter_phi()),5,2))
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
  plt.text(config.N-1,1,'$R_0 < 1$',fontsize=14,color='w',va='top',ha='right')
  if save and config.save:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def gen_1d_plot(data,output,select,save=None,regions=False,legend=None):
  figsize = {'paper': (4,3),                 'isstdr':(4.5,3)}[config.context]
  axespos = {'paper': [0.16,0.14,0.82,0.84], 'isstdr':[0.3,0.16,0.68,0.82]}[config.context]
  plt.figure(figsize=figsize)
  make_1d(data,output,select,regions=regions)
  plt.gca().set_position(axespos)
  if legend:
    plt.legend(legend)
  if save and config.save:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def gen_plots(save=False):
  if config.context == 'paper':
    gen_plots_paper()
  if config.context == 'isstdr':
    gen_plots_isstdr()

def gen_plots_isstdr():
  for output,select in [('prevalence','high'),('prevalence','low')]:
    data = load_data(output,select,tau=[TAU1D])
    gen_1d_plot(data,output,select,fname_fig('isstdr',output,select,tau=TAU1D))

def gen_plots_paper():
  # incidence and prevalence
  for output in OUTPUTS:
    for select in SELECTORS:
      data = load_data(output,select,tau=[TAU1D])
      gen_1d_plot(data,output,select,fname_fig('1d',output,select,tau=TAU1D),regions=True)
      for norm in [False,True]:
        gen_2d_plot(load_data(output,select,norm),fname_fig('surface',output,select,norm=norm))
  # incidence factors
  data = {
    '1d': {output: load_data(output,'I',tau=[TAU1D]) for output in ['X','C'] },
    '2d': {output: load_data(output,'I')             for output in ['X','C'] },
  }
  data['1d']['XC'] = data['1d']['X'] * data['1d']['C'] / CA
  data['2d']['XC'] = data['2d']['X'] * data['2d']['C'] / CA
  for output in ['X','C','XC']:
    gen_1d_plot(data['1d'][output],output,None,fname_fig('1d',output,'I',tau=TAU1D),regions=True)
    gen_1d_plot(data['2d'][output],output,None,fname_fig('2d',output,'I'),regions=True)
  gen_1d_plot(load_data('incidence','all'),'incidence',None,fname_fig('2d','incidence','all'),regions=True)
  # proportion of infections from turnover
  legend = ['High-Risk','Medium-Risk','Low-Risk']
  data = np.concatenate([load_data('tip',select,tau=[TAU1D]) for select in ['high','med','low']])
  gen_1d_plot(data,'tip',None,fname_fig('2d','tip','all',tau=TAU1D),legend=legend)
  

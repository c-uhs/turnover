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

OUTPUTS   = ['prevalence','incidence']
SELECTORS = ['all','high','med','low']
HEALTH    = ['S','I','T']
TAU1D = 0.1
CA = 3 # HACK
# config.save = True

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

def iter_phi_mat(plims=None):
  for phi in iter_phi(plims=plims):
    yield get_sim(phi,tau=TAU1D,tmax=1).params['phi']

# filenames for saving (data, figures)

def fname_fun(*args,**kwargs):
  args = list(args)
  keys = sorted(list(kwargs.keys()))
  for key in keys:
    value = kwargs[key]
    if value is None:
      pass
    elif isinstance(value,bool):
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
  path = os.path.join(config.path['figs'],'sensitivity')
  utils.makedir(path)
  return os.path.join(path,fname_fun(*args,**kwargs)+'.pdf')

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
  p1 = (1/model.params['dur'].iselect(ii=['H'],ki=['M']) - model.params['mu'])/2
  model.params['phi'].update(p1,ii=['H'],ip=['M'])
  model.params['phi'].update(p1,ii=['H'],ip=['L'])
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
  if idx == 'debug':
    sim = run_sim(phi=0.03,tau=TAU1D)
    sim.plot(
        output = 'prevalence',
        selectors = [sim.model.select[s] for s in SELECTORS],
      )
  elif int(idx) == -1:
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
  for select in iter_selectors(sim):
    for output in OUTPUTS:
      save_element(sim.outputs[output],select,tmax,phi,tau)
    for health in HEALTH:
      save_element(sim.outputs['X'],sim.model.select[health].union(select),tmax,phi,tau)
    save_element(sim.outputs['tip'],select,tmax,phi,tau)
  save_element(sim.outputs['C'],sim.model.select['I'],tmax,phi,tau)
  return sim

def save_element(output,select,tmax,phi,tau):
  if config.save:
    fname = fname_element(output.name,select.name,phi=phi,tau=tau)
    value = modelutils.taccum(output,**select).islice(t=tmax)
    utils.savetxt(fname,float(value))

# make surface plots of the result

def load_element(output,select,phi,tau,norm,tol=1e-5):
  phin = list(iter_phi())[0]
  x = np.float(utils.loadtxt(fname_element(output,select,phi=phi,tau=tau)))
  n = np.float(utils.loadtxt(fname_element(output,select,phi=phin,tau=tau))) if norm else 1
  return (x / n) if x > tol else 0

def load_data(output,select,norm=False,phi=None,tau=None,tol=1e-5):
  phi = phi if phi is not None else list(iter_phi())
  tau = tau if tau is not None else list(iter_tau())
  dname = fname_array(output,select,norm=norm)
  if True: # not os.path.exists(dname): # TODO: whats going on here?
    array = [[ load_element(output,select,p,t,norm,tol=tol) for p in phi] for t in tau]
    utils.savecsv(dname,array,append=False)
  else:
    array = utils.loadcsv(dname,asdict=False)
  array = np.array(array,dtype=np.float)
  return array

def make_1d_pretty(output,select=None,health=None):
  ytitle = {
    'prevalence':       'Prevalence',
    'incidence':        'Incidence (per 1000 PY)',
    'X':                'Proportion (\%)',
    'C':                'Average contact rate of infectious',
    'XC':               'Infectious partnerships proportion',
    'tip':              'Proportion of new infections from turnover',
    'dX':               'Rate of change (\%)\n',
    'prevalence-ratio': 'Prevalence ratio:',
    'incidence-ratio':  'Incidence ratio:',
    None: '',
  }
  yspec = {
    'high':     ' among High-Risk',
    'med':      ' among Medium-Risk',
    'low':      ' among Low-Risk',
    'all':      ' Overall',
    'high-low': ' High vs Low Risk',
    'high-med': ' High vs Medium Risk',
    'med-low':  ' Medium vs Low Risk',
    None: '',
  }
  yhspec = {
    'S': ' Susceptible',
    'I': ' Infectious',
    'T': ' Treated',
    None: '',
  }
  ylabel = ytitle[output]+yspec[select]+yhspec[health]
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

def make_1d(data,output,select,health=None,regions=False,cmap=None):
  def get_peaks(data,interp):
    x  = np.arange(0,config.N,1)
    xi = np.arange(0,config.N,interp)
    fi = lambda y: interpolate.interp1d(x,y,kind='cubic',bounds_error=False)(xi)
    xp = [np.nanargmax(fi(y))*interp for y in data]
    yp = [np.nanmax(fi(y)) for y in data]
    return xp,yp
  if cmap is None:
    cmap = plt.cm.Blues_r(np.linspace(0,0.8,data.shape[0]))
  plt.gca().set_prop_cycle('color',cmap)
  plt.plot(range(0,config.N),data.transpose())
  if regions:
    xp,yp = get_peaks(data,0.01)
    if np.any(xp):
      draw_regions(xp,yp,['A','B'])
  make_1d_pretty(output,select,health)

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

def gen_1d_plot(data,output,select,health=None,save=None,regions=False,legend=None,cmap=None):
  figsize = {'paper': (4,3),                 'isstdr':(4.5,3)}[config.context]
  axespos = {'paper': [0.17,0.14,0.81,0.84], 'isstdr':[0.3,0.16,0.68,0.82]}[config.context]
  plt.figure(figsize=figsize)
  make_1d(data,output,select,health=health,regions=regions,cmap=cmap)
  plt.gca().set_position(axespos)
  if legend:
    plt.legend(legend)
  if save and config.save:
    plt.savefig(save)
  else:
    plt.show()
  plt.close()

def gen_flows_plot(health,select):
  sim = get_sim(phi=0.1,tau=TAU1D,tmax=1)
  params = sim.params
  phix = np.stack(iter_phi_mat())
  sall = ['high','med','low']
  # health state flows
  px  = params['px'].iselect(**sim.model.select[select])
  pxh   = load_data('X',health+' '+select,tau=[TAU1D])
  inc   = load_data('X','S '+select,tau=[TAU1D]) * load_data('incidence',select,tau=[TAU1D]) / 1000
  treat = load_data('X','I '+select,tau=[TAU1D]) * TAU1D
  death = pxh * params['mu']
  birth = px  * params['nu'] * np.ones(death.shape)
  # turnover flows
  ii = sall.index(select)
  io = [i for i in [0,1,2] if i is not ii]
  eff = load_data('X',health+' '+select,tau=[TAU1D]) * np.sum(phix[:,ii,:],axis=1)
  aff = np.sum([
      load_data('X',health+' '+sall[i],tau=[TAU1D]) * phix[:,i,ii]
      for i in io ], axis=0)
  specs = [
      ('Birth',         birth, [0.0,0.6,1.0]),
      ('Incidence',     inc,   [0.8,0.2,0.0]),
      ('Treatment',     treat, [0.8,0.2,0.6]),
      ('Death',         death, [0.4,0.4,0.4]),
      ('Turnover into', aff,   [1.0,0.6,0.0]),
      ('Turnover from', eff,   [1.0,0.8,0.0]),
    ]
  if health == 'S':
    specs.pop(2)
    signs = [+1,-1,-1,+1,-1]
  if health == 'I':
    specs.pop(0)
    signs = [+1,-1,-1,+1,-1]
  if health == 'T':
    specs.pop(0)
    specs.pop(0)
    signs = [+1,-1,+1,-1]
  data = np.concatenate(tuple(spec[1]*sign for spec,sign in zip(specs,signs))) / pxh * 100
  # data = np.concatenate((data,np.sum(data,axis=0,keepdims=True)))
  gen_1d_plot(data,'dX',select,health,
      save   = fname_fig('dX',select,health),
      legend = [spec[0] for spec in specs],
      cmap   = [spec[2] for spec in specs],
    )

def gen_plots(save=False):
  if config.context == 'paper':
    gen_plots_paper()
  if config.context == 'isstdr':
    gen_plots_isstdr()

def gen_plots_isstdr():
  for output,select in [('prevalence','high'),('prevalence','low')]:
    data = load_data(output,select,tau=[TAU1D])
    gen_1d_plot(data,output,select,save=fname_fig('isstdr',output,select,tau=TAU1D))

def gen_plots_paper():
  pass
  # incidence and prevalence
  # for output in OUTPUTS:
  #   for select in SELECTORS:
  #     data = load_data(output,select,tau=[TAU1D])
  #     gen_1d_plot(data,output,select,save=fname_fig('1d',output,select,tau=TAU1D),regions=True)
  #     for norm in [False,True]:
  #       gen_2d_plot(load_data(output,select,norm),save=fname_fig('surface',output,select,norm=norm))
  #   for pair in [('high','med'),('high','low'),('med','low')]:
  #     select = '-'.join(pair)
  #     data = load_data(output,pair[0],tau=[TAU1D],tol=0) / load_data(output,pair[1],tau=[TAU1D],tol=0)
  #     gen_1d_plot(data,output+'-ratio',select,save=fname_fig('1d','ratio',output,select,tau=TAU1D))
  #     data = load_data(output,pair[0],tol=0) / load_data(output,pair[1],tol=0)
  #     gen_1d_plot(data,output+'-ratio',select,save=fname_fig('2d','ratio',output,select))
  # # health states
  # cmap = [[0.0,0.6,1.0],[0.8,0.2,0.0],[0.8,0.2,0.6]]
  # legend = ['Susceptible','Infectious','Treated']
  # for select in ['all']:#SELECTORS:
  #   data = np.concatenate([
  #       load_data('X',health+' '+select,tau=[TAU1D])
  #       for health in HEALTH
  #     ])
  #   # data /= np.sum(data,axis=0)
  #   data -= np.mean(data,axis=1).reshape((3,1))
  #   gen_1d_plot(100*data,'X',select,None,save=fname_fig('1d','X','health',select,tau=TAU1D),cmap=cmap,legend=legend)
  # # incidence factors
  # ositer = [('C','I'),('prevalence','all'),('incidence','all')]
  # dtiter = [('1d',[TAU1D]),('2d',None)]
  # for dim,tau in dtiter:
  #   data = { output: load_data(output,select,tau=tau) for output,select in ositer }
  #   data['XC'] = data['prevalence'] * data['C'] / CA
  #   for output,select in ositer+[('XC','I')]:
  #     gen_1d_plot(data[output],output,None,save=fname_fig(dim,output,select,tau=tau),regions=True)
  # # overall incidence
  # gen_1d_plot(load_data('incidence','all'),'incidence',None,save=fname_fig('2d','incidence','all'),regions=True)
  # equilibrium flows
  # for health in HEALTH:
  #   for select in ['high','med','low']:
  #     gen_flows_plot(health,select)
  # # proportion of infections from turnover
  # legend = ['High-Risk','Medium-Risk','Low-Risk']
  # data = np.concatenate([load_data('tip',select,tau=[TAU1D]) for select in ['high','med','low']])
  # gen_1d_plot(data,'tip',None,fname_fig('2d','tip','all',tau=TAU1D),legend=legend)
  

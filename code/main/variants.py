import os,sys
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
config.epimodel()
config.plot()
import numpy as np
config.numpy()
import matplotlib.pyplot as plt
from collections import OrderedDict as odict
from textwrap import wrap
from copy import deepcopy
import utils
import modelutils
import system

CONTEXT = 'isstdr' # in ['paper','isstdr']
SAVE = True

def fname_fig(compare,output,selector,**params):
  if CONTEXT == 'paper':
    path = [config.path['figs'],'plots','compare']
  if CONTEXT == 'isstdr':
    path = [config.path['root'],'docs','conferences','isstdr-2019','slides','figs']
  return os.path.join(*path,
    '-'.join(
      [compare,output,selector]+
      ['{}={}'.format(name,value) for name,value in params.items()])+'.pdf'
  )

def load_fit(name,sim):
  sim = deepcopy(sim)
  model = sim.model
  fname = os.path.join(config.path['data'],'fit',shortname(name)+'.json')
  model.params.fromdict(utils.loadjson(fname))
  sim.init_model(model)
  sim.init_params()
  return sim

def get_sim(variant=None,t=None):
  specs = system.get_specs()
  model = system.get_model()
  if variant in [None,0,'full']:
    pass
  elif variant in [1,'no-hetero']:
    model.collapse(['ii'])
  elif variant in [2,'no-growth']:
    model.params['nu'].update(model.params['mu'])
  elif variant in [3,'no-turnover']:
    model.params['dur'].update(np.nan)
    model.params['phi'].update(np.nan)
  return system.get_simulation(model,t=t)

def shortname(name):
  short = name.split(' ')[0]
  return short if 'fit' not in name else short+'-fit'

def txtsave(name,sim,output,selector,txt):
  fname = os.path.join(config.path['data'],'fit','-'.join([name,output,selector])+'.txt')
  value = modelutils.taccum(sim.outputs[output],**sim.model.select[selector]).islice(t=sim.t[-1])
  utils.savetxt(fname,txt(value) if callable(txt) else str(value))

def plot_iter(sims,output,selector,txt=False):
  legend = []
  colors = [[0.8,0.0,0.0],[1.0,0.6,0.6],[0.8,0.0,0.0],[1.0,0.6,0.6]]
  linestyles = ['-','-','--','--']
  ylim = {'paper': None, 'isstdr': [0.5, 1.0]}[CONTEXT]
  for (name,sim),color,ls in zip(sims.items(),colors,linestyles):
    legend.append(name)
    select = sim.model.select[selector]
    select.color = color
    sim.plot(
      output = output,
      selectors = [sim.model.select[selector]],
      xlabel = 'Time (years)',
      show = False,
      legloc = False,
      linestyle = ls,
      ylim = ylim,
    )
    if txt:
      txtsave(shortname(name),sim,output,selector,txt)
  if CONTEXT == 'paper':
    plt.legend(legend)
  if CONTEXT == 'isstdr':
    plt.legend(['Turnover','No Turnover'], loc='lower right')
    plt.xlabel(plt.gca().get_xlabel(),fontsize='x-large')
    plt.ylabel('\n'.join(wrap(plt.gca().get_ylabel(),13,break_on_hyphens=False)),
      fontsize='x-large',
      labelpad=40,
      rotation=0,
      va='center',
    )

def run_sim(sim,outputs=None):
  outputs = outputs if outputs is not None else []
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = outputs
  ))
  return sim.solve()

def exp_run_plot(compare,sims,outputs,selectors,save=False,txt=False,**params):
  for sim in sims.values():
    for name,value in params.items():
      if name in sim.model.params:
        sim.model.params[name].update(value)
    sim.update_params(sim.model.params)
    run_sim(sim,outputs)
  figsize = {'paper': (4,3), 'isstdr': (4.5,3)}[CONTEXT]
  axespos = {'paper': [0.16,0.14,0.82,0.84], 'isstdr':[0.33,0.16,0.65,0.82]}[CONTEXT]
  for output in outputs:
    for selector in selectors:
      plt.figure(figsize=figsize)
      plot_iter(sims,output,selector,txt=txt)
      plt.gca().set_position(axespos)
      if save and SAVE:
        plt.savefig(fname_fig(compare,output,selector,**params))
        plt.close()
      else:
        plt.show()

def exp_hetero(save=False):
  sims = odict([
    ('Full (Risk Heterogeneity)',  get_sim('full')),
    ('V1 (No Risk Heterogeneity)', get_sim('no-hetero')),
  ])
  exp_run_plot('hetero',
    sims      = sims,
    outputs   = ['prevalence','incidence'],
    selectors = ['all'],
    save      = save,
  )

def exp_growth(save=False):
  sims = odict([
    ('Full (Population Growth)',  get_sim('full')),
    ('V2 (No Population Growth)', get_sim('no-growth')),
  ])
  exp_run_plot('growth',
    sims      = sims,
    outputs   = ['prevalence','incidence'],
    selectors = ['all','high','low'],
    save      = save,
  )

def exp_turnover(save=False):
  sims = odict([
  ('Full (Turnover)',  get_sim('full')),
  ('V3 (No Turnover)', get_sim('no-turnover')),
])
  for tau in [0.1, 0.2]:
    exp_run_plot('turnover',
      sims      = sims,
      outputs   = ['prevalence','incidence'],
      selectors = ['all','high','med','low'],
      save      = save,
      tau       = tau,
      # infect    = [[0.5/3],[2.0/3],[7.5/3]], # TEMP
    )

def exp_tpaf(save=False):
  tmax = {'paper': 100, 'isstdr': 30}[CONTEXT]
  t = system.get_t(tmax=tmax)
  for case in ['raw','fit','both']:
    print(case,flush=True)
    sims = odict([
      ('Full (Turnover)',  get_sim('full',t=t)),
      ('V3 (No Turnover)', get_sim('no-turnover',t=t)),
    ])
    if case == 'raw':
      pass
    if case == 'fit':
      for name in list(sims.keys()):
        sims.update([(name, load_fit(name,sims[name]))])
    if case == 'both':
      for name in list(sims.keys()):
        sims.update([(name+' [fit]', load_fit(name,sims[name]))])
    exp_run_plot('tpaf',
      sims      = sims,
      outputs   = ['tpaf-high'],
      selectors = ['all'],
      save      = save,
      vs        = case,
    )
  # equilibrium prevalence plot
  # names = list(sims.keys())
  # for name in names:
  #   sim_eq = sims[name].model.equilibriate(tmax=500,tol=1e-6)
  #   sims[name]._model.X0 = sim_eq.X.islice(t=sim_eq.teq)
  #   sims[name].model.params['infect'].update(0)
  # exp_run_plot('tpaf',
  #   sims      = sims,
  #   outputs   = ['prevalence'],
  #   selectors = ['low','high'],
  #   save      = save,
  #   txt       = lambda x: '{:.0f}\%'.format(100*float(x)),
  # )

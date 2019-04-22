import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot(tex=True)

import numpy as np
config.numpy()
from collections import OrderedDict as odict
import matplotlib.pyplot as plt
import elements

import system
import variants
import utils
import batch

def fname_fig(name):
  return os.path.join(config.path['figs'],'plots','compare',name)
  
def load_fit(name,sim):
  model = sim.model
  model.params.fromdict(utils.loadjson(batch.fname_fit(name)))
  sim.init_model(model)
  sim.init_params()
  return sim

def add_plot(sim,selector,output):
  # initialize outputs
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = [output]
  ))
  ylabel = '{}{}{}'.format(
    'TPAF of High Activity Women' if output=='tpaf-WH' else output.capitalize(),
    ' among {}'.format(selector.title) if selector else '',
    ' (per 1000)' if output == 'incidence' else ''
  )
  sim.solve()
  sim.plot(
    output = output,
    selectors = [selector],
    xlabel = 'Time (years)',
    # ylabel = ylabel, # TEMP
    show = False,
    leg = False,
  )

def make_plot(output,sim_iter,save=None,t=None,legloc=None):
  legend = []
  for title,sim,selector in sim_iter(t=t):
    add_plot(sim,selector,output)
    legend += [title]
  legloc = legloc if legloc is not None else 'lower right'
  plt.legend(legend,loc=legloc,fontsize=8)
  if save:
    plt.savefig(fname_fig(save))
    plt.close()
  else:
    plt.show()

def select_iter():
  for color in [[0.8,  0,  0],[  0,  0,0.8]]:
    yield elements.Selector(
    name = None,
    select = {},
    title = None,
    color = color,
    linestyle = '-',
  )

def growth():
  def sim_iter(t=None):
    selecti = select_iter()
    for i,(name,sim) in enumerate(variants.get_sims(t=t).items()):
      if i in [0,3]:
        selector = next(selecti)
        select = sim.model.select['all']
        selector.update(select)
        title = variants.make_title(*variants.parse_name(name),i=i)
        yield title,sim,selector
  make_plot(
    output = 'prevalence',
    sim_iter = sim_iter,
    save = 'compare-growth-prevalence.eps',
  )

def hetero():
  def sim_iter(t=None):
    selecti = select_iter()
    for i,(name,sim) in enumerate(variants.get_sims(t=t).items()):
      if i in [1,2]:
        selector = next(selecti)
        title = variants.make_title(*variants.parse_name(name),i=i)
        yield title,sim,selector
  make_plot(
    output = 'prevalence',
    sim_iter = sim_iter,
    save = 'compare-hetero-prevalence.eps',
  )

def turnover(*outputs):
  def sim_iter(who,t=None):
    selecti = select_iter()
    for i,(name,sim) in enumerate(variants.get_sims(t=t).items()):
      if i in [0,1]:
        selector = next(selecti)
        select = sim.model.select[who]
        selector.update(select)
        title = variants.make_title(*variants.parse_name(name),i=i)+' '+select.title
        yield title,sim,selector
  for output in outputs:
    for who in ['all','high','low']:
      make_plot(
        output = output,
        sim_iter = lambda t: sim_iter(who=who,t=t),
        save = 'compare-turnover-{}-{}.eps'.format(output,who),
      )

def tpaf():
  def sim_iter(t=None):
    t = system.get_t(tmax=100)
    selecti = select_iter()
    for i,(name,sim) in enumerate(variants.get_sims(t=t).items()):
      if i in [0,1]:
        selector = next(selecti)
        title = variants.make_title(*variants.parse_name(name),i=i)
        yield title,sim,selector
        sim = load_fit(name,sim)
        title += ' [fit]'
        selector.specs.update(linestyle = '--')
        yield title,sim,selector
  make_plot(
    output = 'tpaf-WH',
    sim_iter = sim_iter,
    save = 'compare-turnover-tpaf-fit.eps',
  )

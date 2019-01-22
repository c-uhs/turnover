import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot(tex=True)

import numpy as np
config.numpy()
from collections import OrderedDict as odict
import matplotlib.pyplot as plt

import system
import variants
import modelutils
import utils

def savename(*args):
  return os.path.join(config.path['figs'],'plots',*args)

def loadfit(sim,fname):
  model = sim.model
  model.params.fromdict(utils.loadjson(fname))
  sim.init_model(model)
  sim.init_params()
  return sim

def iter_tpaf(fitdir):
  t = system.get_t(tmax=100)
  for i,(name,sim,selector) in enumerate(variants.iter_all(t=t)):
    nu,G,Z = variants.parse_name(name)
    if (G != 1) and (nu != sim.model.params['mu']):
      selector.title = 'V{}: {}'.format(i+1,selector.title)
      yield name,sim,selector
  for i,(name,sim,selector) in enumerate(variants.iter_all(t=t)):
    nu,G,Z = variants.parse_name(name)
    if (G != 1) and (nu != sim.model.params['mu']):
      sim = loadfit(sim,os.path.join(fitdir,name+'.json'))
      selector.specs['linestyle'] = '--'
      selector.title = 'V{}: {} [fit]'.format(i+1,selector.title)
      yield name,sim,selector

if __name__ == '__main__':
  select = {}
  legend = []
  output = 'tpaf-WH'
  fitdir = os.path.join(config.path['specs'],'fit','C')
  for name,sim,selector in iter_tpaf(fitdir):
    selector.update(select)
    sim.init_outputs(system.get_outputs(
      spaces = sim.model.spaces,
      select = sim.model.select,
      t = sim.t,
      names = [output]
    ))
    print(name,flush=True)
    sim.solve()
    sim.plot(
      outputs = [output],
      selectors = [selector],
      xlabel = 'Time (years)',
      show = False,
      leg = False,
      ylabel = 'TPAF of High Risk Women',
    )
    legend += [selector.title]
  plt.legend(legend,loc='lower right')
  plt.savefig(savename('compare','tpaf-fit-not-C-compare.eps'))
  plt.show()
  plt.close()
  
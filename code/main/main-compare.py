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

def make_plot(output,save,fitdir=None,t=None,select={},**kwargs):
  legend = []
  for i,(name,sim,selector) in enumerate(variants.iter_all(t=t)):
    if fitdir:
      fname = os.path.join(fitdir,name+'.json')
      if not os.path.exists(fname): continue
      model = sim.model
      model.params.fromdict(utils.loadjson(fname))
      sim.init_model(model)
      sim.init_params()
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
      **kwargs,
    )
    legend += ['V{}: {}'.format(i+1,selector.title)]
  plt.legend(legend,loc='lower right')
  plt.savefig(save)
  plt.show()
  plt.close()

if __name__ == '__main__':
  # make_plot(output = 'prevalence') # testing
  make_plot(
    output = 'tpaf-WH',
    save = savename('compare','tpaf-C-compare.eps'),
    # fitdir = os.path.join(config.path['specs'],'fit','C'),
    t = system.get_t(tmax = 200),
    ylabel = 'TPAF of High Risk Women',
  )
  # selectors = system.get_specs()['select']
  # for name in ['all','high','med','low']:
  # # for name in ['S','I','T']:
  #   print(name,flush=True)
  #   make_plot(
  #     output = 'prevalence',
  #     save = savename('compare','prevalence-{}-compare.eps'.format(name)),
  #     select = selectors[name],
  #     t = system.get_t(tmax = 200),
  #     ylabel = 'Prevalence among {}'.format(selectors[name].title),
  #   )
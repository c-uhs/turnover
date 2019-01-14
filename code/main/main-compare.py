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

def main_plot(output,**kwargs):
  legend = []
  for i,(name,sim,selector) in enumerate(variants.iter_all()):
    sim.init_outputs(system.get_outputs(
      spaces = sim.model.spaces,
      select = sim.model.select,
      t = sim.t,
      names = [output]
    ))
    sim.solve()
    sim.plot(
      outputs = [output],
      selectors = [selector],
      xlabel = 'Time (years)',
      show = False,
      leg = False,
      **kwargs,
    )
    legend += ['V{}: '.format(i+1)+selector.title]
  plt.legend(legend)
  plt.savefig(savename('compare',output+'-compare.eps'))
  plt.close()

if __name__ == '__main__':
  main_plot(
    output = 'incidence',
    ylabel = 'Overall Incidence (per 1000)',
  )
  main_plot(
    output = 'prevalence',
    ylabel = 'Overall Prevalence',
  )

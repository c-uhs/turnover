# load code / config
import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
config.epimodel()
config.numpy()
config.plot()
# external module
import re
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
# epi-model modules
import utils
import modelutils
from space import Space,Array
from elements import Color
# relative imports
import system
import sensitivity

names = ['S','I','R']
# names = ['S','infected']
lights = [0.6,0.0,0.3]
out = 'X'
n = 1

def make_pie(sim,select,phi):
  plt.figure(figsize=(2,2))
  selectors = [
    sim.model.select[name].union(sim.model.select[select])
    for name in names
  ]
  SIR = np.array([
    modelutils.taccum(sim.outputs[out],**selector).islice(t=sim.t[-1])
    for selector in selectors
  ])
  SIR = SIR / SIR.sum()
  colors = [selector.color.lighten(light) for selector,light in zip(selectors,lights)]
  labels = [sim.model.select[name].label for name in names]
  reorder = lambda x: [x[0],x[2],x[1]]
  # reorder = lambda x: [x[0],x[1]]
  plt.pie(reorder(SIR), colors=reorder(colors), startangle=90, counterclock=True )
  plt.tight_layout(pad=-1.8)
  figdir = os.path.join(config.path['figs'],'flows','phi={}'.format(phi))
  utils.makedir(figdir)
  if config.save:
    plt.savefig(os.path.join(figdir,'{}-{}.pdf'.format('flow',select,phi)),transparent=True)
  else:
    plt.show()
  plt.close()
  make_legend(labels,colors)

def make_legend(labels,colors):
  # plt.figure(figsize=(2.2,0.5))
  plt.figure(figsize=(1.5,1))
  for color in colors:
    plt.fill([np.nan]*3,[np.nan]*3,color=color)
  # plt.legend(['$\\mathcal{S}$','$\\mathcal{I}$','$\\mathcal{R}$'],loc='center',ncol=3)
  plt.legend(labels,loc='center')
  plt.gca().set_axis_off()
  if config.save:
    plt.savefig(os.path.join(config.path['figs'],'flows','flows-legend.pdf'))
  else:
    plt.show()
  plt.close()

def make_tikz(label,phi):
  tikzdir = os.path.join(config.path['tikz'],'flows');
  phistr  = 'phi={}'.format(phi)
  flowdir = os.path.join(config.path['figs'],'flows')
  # What is this 3x escape mess?
  configstr = \
    '\\\\graphicspath{{'+os.path.join(flowdir,phistr)+'/}}'+\
    '\\\\\\\\newcommand{\\\\\\\\turnover}{'+str(4*np.minimum(1,phi**(1/3)))+'}'
  os.system('cd {} && echo {} > config.tex && pdflatex flows.tex >/dev/null && cp flows.pdf {}/{}'.format(
    tikzdir, configstr, flowdir, 'flow-{}.pdf'.format(label) ))

def make_figs():
  phis = list(sensitivity.iter_phi())
  for label,phi in [
      ('low',    phis[n]),
      ('med',    phis[int((config.N-1)/2)]),
      ('high',   phis[config.N-n-1]),
      ('extreme',10),
    ]:
    specs = system.get_specs()
    model = system.get_model()
    sim = sensitivity.get_sim(phi,0.1)
    sim.init_outputs(system.get_outputs(
      spaces = sim.model.spaces,
      select = sim.model.select,
      t = sim.t,
      names = [out]
    ))
    if label == 'extreme':
      sim.update_params(dict(ibeta=0.038))
    sim.solve()
    for name in ['high','low']:
      make_pie(sim,name,phi)
    make_tikz(label,phi)

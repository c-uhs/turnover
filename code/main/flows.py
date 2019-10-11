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

names = ['S','I','T']
selects = ['high','med','low']
# names = ['S','infected']
lights = [0.0,0.0,0.0]
out = 'X'
n = 1

def gen_pie_data(sim,select,phi):
  SIR = np.array([
    modelutils.taccum(
      sim.outputs[out],
      **sim.model.select[name].union(sim.model.select[select])
    ).islice(t=sim.t[-1])
    for name in names
  ])
  SIR = np.round(SIR * 360 / SIR.sum()).astype(np.int)
  SIR[0] = 360-SIR[1:].sum() # HACK in case sum < 1
  tdir = os.path.join(config.path['data'],'flows','phi={}'.format(phi))
  if config.save:
    utils.makedir(tdir)
    for name,value in zip(names,SIR):
      utils.savetxt(os.path.join(tdir,'flow-{}-{}.tex'.format(select,name)),int(np.round(value)))

def make_tikz(label,phi):
  tikzdir = os.path.join(config.path['tikz'],'flows')
  tdir    = os.path.join(config.path['data'],'flows','phi={}'.format(phi))
  flowdir = os.path.join(config.path['figs'],'flows')
  # What is this 3x escape mess?
  configstr = '\n'.join([
    '\\newcommand{{\\x{}{}}}{{{}}}'.format(name,select,
        int(utils.loadtxt(os.path.join(tdir,'flow-{}-{}.tex'.format(select,name))))
      ) for name in names for select in selects
  ])+'\n\\newcommand{\\turnover}{'+str(4*np.minimum(1,phi**(1/3)))+'}'
  utils.savetxt(os.path.join(tikzdir,'config.tex'),configstr)
  os.system('cd {} && pdflatex flows.tex >/dev/null && cp flows.pdf {}/{}'.format(
    tikzdir, flowdir, 'flows-{}.pdf'.format(label) ))

def run_sims():
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
    for select in selects:
      gen_pie_data(sim,select,phi)
    make_tikz(label,phi)

# load code / config
import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
config.epimodel()
config.numpy()
config.plot()
# external module
import numpy as np
import matplotlib.pyplot as plt
# epi-model modules
from utils import relimport
from modelutils import partner
# relative imports
import variants
import system
import surface

N = 9

def main():
  # sims = variants.get_sims(t=system.get_t(tmax=100))
  # sim = sims['nu=0.05_G=3_Z=1']
  out = 'inf-ratio'
  colors = plt.get_cmap('plasma',N).colors
  for dix in [5,10,20,50]:
    legend = []
    for idh,idi,dh,di in surface.iter_both(N):
      if idi == 0:
        print(dh,flush=True)
        legend += ['$\delta_H = {}$'.format(dh)]
        sim = surface.get_sim(dh,dix)
        sim.init_outputs(system.get_outputs(
          spaces = sim.model.spaces,
          select = sim.model.select,
          t = sim.t,
          names = [out]))
        sim.solve()
        selector = sim.model.select['WL']
        selector.update(ip=['M','H'])
        selector.color = colors[idh]
        sim.plot(
          output = out,
          selectors = [selector],
          show = False,
          leg = False,
        )
    plt.ylim([0,1])
    plt.legend(legend)
    plt.ylabel('Proportion of New Infections among Low Activity from Turnover')
    plt.savefig('prop-IT-in-L-from-HM-di={:02d}.pdf'.format(dix))
    plt.close()
  # plt.show()

  # sim.init_outputs(system.get_outputs(
  #   spaces = sim.model.spaces,
  #   select = sim.model.select,
  #   t = sim.t,
  #   names = [out]))
  # sim.solve()
  # selectors = [sim.model.select[name] for name in ['WH','WM','WL']]
  # for selector in selectors:
  #   selector.update({'ip':selector.pop('ii')})
  #   selector.update({'ii':'L'})
  # sim.plot(
  #   output = out,
  #   selectors = selectors,
  # )

import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
import numpy as np
config.epimodel()
config.plot(tex=True)
import matplotlib.pyplot as plt
from collections import OrderedDict as odict
import modelutils
import elements
import system

output = 'prevalence'
N = 7

def idur(i):
  return 10**-i

def iterdur():
  for i in np.logspace(-1.5,+1.5,N)[::-1]:
    yield i

def get_sim(type,dh=None):
  model = system.get_model()
  if type == 'homo':
    model.collapse(['ii','ip'])
    name = '$G = 1$ (Homogeneous)'
    dt = 0.2
  if type == 'zeta':
    dhmax = 1.0/model.params['mu']
    model.params['pe'].update(np.nan)
    model.params['zeta'].update(np.nan)
    model.params['dur'].update([dh,min(dh*5,dhmax),min(dh*30,dhmax)])
    name = '$G = 3, \\delta_H = {:.03f}$'.format(model.params['dur'][0,0])
    dt = min(dh*0.8,0.2)
  t = system.get_t(dt=dt,tmin=0,tmax=200)
  sim = system.get_simulation(model,t=t)
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = [output],
  ))
  return sim,name

def plotfun(sim,name,color,**specs):
  selector = sim.model.select['all']
  selector.color = color
  selector.specs.update(**specs)
  sim.solve()
  sim.plot(
    outputs = [output],
    selectors = [selector],
    ylabel = 'Prevalence Overall',
    show = False,
    leg = False,
  )
  return [name]

if __name__ == '__main__':
  legend = []
  cmap = elements.Color([1,0,0]).cmap(N,[-0.7 ,+0.7])[::-1]
  for dh,c in zip(iterdur(),cmap):
    print(dh,flush=True)
    legend += plotfun(*get_sim('zeta',dh=dh),color=c)
  legend += plotfun(*get_sim('homo'),color=[0,0,0],linestyle='--')
  plt.legend(legend,loc='lower right',fontsize=8)
  plt.savefig(os.path.join(config.path['figs'],'plots','compare','homo-vs-zeta-prev-full.eps'))
  plt.show()


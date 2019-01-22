import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
import numpy as np
config.epimodel()
config.plot(tex=True)
import matplotlib.pyplot as plt
from collections import OrderedDict as odict
import utils
import run2d

def make_plot(output,select,N,norm=False):
  def clean(data,norm=False):
    return [[
      np.float(d) / (1 if not norm else np.float(da[0]))
      for d in da]
      for da in data]
  def ticks(values,incr,dec):
    return range(len(values))[::incr],[np.around(v,dec) for v in values][::incr]
  def fname_fig(output,select,norm):
    return os.path.join(config.path['figs'],'2d',
      '2d-{}{}-{}.eps'.format(output,'-norm' if norm else '',select))
  data = clean(utils.loadcsv(run2d.fname_data(output,select,N),asdict=False),norm)
  plt.figure(figsize=(6.5,5))
  plt.imshow(data,cmap=plt.get_cmap('inferno'),interpolation='none')
  plt.colorbar().ax.tick_params(labelsize=16)
  plt.yticks(*ticks(list(run2d.iter_tau(N)),2,2),fontsize=16)
  plt.xticks(*ticks(list(run2d.iter_dh(N)),2,1),fontsize=16)
  plt.ylabel('$\\tau$',fontsize=20)
  plt.xlabel('$\\delta_H$',fontsize=20)
  plt.grid(None)
  plt.tight_layout(pad=0.5)
  plt.savefig(fname_fig(output,select,norm))
  # plt.title('{}: {}'.format(output.capitalize(),select.capitalize()))
  # plt.show()
  plt.close()

if __name__ == '__main__':
  N = 15
  for output in ['incidence','prevalence']:
    for select in ['all','high','med','low']:
      make_plot(output,select,N,norm=False)
      make_plot(output,select,N,norm=True)


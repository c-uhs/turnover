import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot()

import numpy as np
config.numpy()
from collections import OrderedDict as odict
from matplotlib.pyplot import close as plt_close

import system
import variants
import modelutils
import calibration
import utils

def runsim(sim,plots,varname,vset):
  # run a single simulation for a single model variant, and generate named plots
  specs = system.get_specs()

  def savename(vset,plot,varname,ftype='png'):
    # filename for the plot
    fname = plot+'_'+varname+'.'+ftype
    return os.path.join(config.path['root'],'outputs','figs','plots',vset,fname)

  def specfun(**spec):
    # replace names with selectors
    spec['selectors'] = [sim.model.select[name] for name in spec['selectors']]
    # assign default plot fun if none assigned
    spec['fun'] = spec.pop('fun',lambda **kwargs: sim.plot(**kwargs))
    return spec

  def tpafplotfun(**spec):
    # special plotting function for tPAF since it is not a model output
    spec.pop('outputs')
    pop = spec.pop('pop')
    tpaf = modelutils.tpaf(sim,pop,dxinf='xlam',beta='beta',copy=True)
    modelutils.plot(sim.t,tpaf,**spec)

  def fitfun(**spec):
    # special case: fit model to targets, save fitted parameters (from params.json)
    # initialize the targets and calsim
    targets = system.get_targets(
      sim.model.spaces,
      sim.model.select,
      sim.outputs,
      t = sim.t)
    calsim = calibration.CalibrationSim('calsim',
      sim = sim,
      targets = targets)
    # initialize the results file if it does not already exist
    fname = os.path.join(config.path['root'],'outputs','params','optimized.json')
    if not os.path.exists(fname):
      with open(fname,'w') as f:
        f.write('{}')
    results = utils.loadjson(fname)
    # run the optimization
    opt = calsim.optimize(plot='tmp.png')
    # update and order the results file
    results.update({spec['title']:{name:param.tolist() for name,param in opt.items()}})
    results = odict([(p,v) for p,v in sorted(results.items())])
    utils.savejson(fname,results,indent=2)

  # dictionary of specifications for the plots
  specs = { name:spec for name,spec in \
    { 'N': specfun(
        outputs = ['N'],
        selectors = ['WH','MH','WM','MM','WL','ML'],
        ylim = [0,1500]),
      'X-sit': specfun(
        outputs = ['X'],
        selectors = ['S','I','T'],
        ylim = [0,1.0]),
      'X-groups': specfun(
        outputs = ['X'],
        selectors = ['WH','MH','WM','MM','WL','ML'],
        ylim = [0,0.55]),
      'prevalence': specfun(
        outputs = ['prevalence'],
        selectors = ['WH','WM','WL'],
        ylim = [0,0.8]),
      'incidence': specfun(
        outputs = ['incidence'],
        selectors = ['WH','WM','WL'],
        ylim = [0,300]),
      'incidence-abs': specfun(
        outputs = ['incidence-abs'],
        selectors = ['WH','WM','WL'],
        ylim = [0,50]),
      'cum-infect': specfun(
        outputs = ['cum-infect'],
        selectors = ['WH','WM','WL'],
        ylim = [0,1000]),
      'tpaf-fsw': specfun(
        fun = tpafplotfun,
        outputs = ['tPAF'],
        pop = {'kp':'W','ip':'H'},
        selectors = ['all'],
        ylabel = 'tPAF of High Risk Group',
        ylim = [0,1]),
      'fit': specfun(
        fun = fitfun,
        outputs = ['prevalence'],
        selectors = [],
      ),
    }.items() if (name in plots) or ((len(plots)==0) and (name != 'fit')) }

  # initialize the required outputs
  outputs = utils.unique(utils.flatten(spec['outputs'] for spec in specs.values()))
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = outputs))
  # solve the system
  print('> {}'.format(varname),flush=True)
  sim.solve()
  # generate and save the specified plots
  for plot,spec in specs.items():
    print('  + {}'.format(plot),flush=True)
    spec.pop('fun')(**spec,
      show = False,
      title = varname,
      save = savename(vset,plot,varname))
    plt_close()

def print_turnover(name,sim):
  print('-'*50)
  print(name)
  print(sim.model.params['zeta'].islice(t=2000,ki='M'))
  print(sim.model.params['dur'].islice(t=2000,ki='M'))
  print(sim.model.params['pe'].islice(t=2000,ki='M'))

if __name__ == '__main__':

  # variants w.r.t. model structure
  for name,sim in variants.get_sims().items():
    # print_turnover(name,sim)
    runsim(sim,[],name,'structure')


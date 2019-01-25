import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot(tex=True)

import numpy as np
config.numpy()
from collections import OrderedDict as odict
from matplotlib.pyplot import close as plt_close

import system
import variants
import modelutils
import calibration
import utils

def fname_fig(plot,varname,ftype='png'):
  # filename for the plot
  fname = plot+'_'+varname+'.'+ftype
  return os.path.join(config.path['figs'],'plots','batch',fname)

def fname_fit(varname):
  return os.path.join(config.path['data'],'fit',varname+'.json')

def run_sim(sim,todo,varname):
  # run a single simulation for a single model variant, and generate named plots
  specs = system.get_specs()

  def specfun(**spec):
    # replace names with selectors
    spec['selectors'] = [sim.model.select[name] for name in spec['selectors']]
    # assign default plot fun if none assigned
    spec['fun'] = spec.pop('fun',lambda **kwargs: sim.plot(**kwargs))
    return spec

  def fitfun(**spec):
    # special case: fit model to targets, save fitted parameters (from params.json)
    # initialize the targets and calsim
    targets = system.get_targets(
      sim.model.spaces,
      sim.model.select,
      sim.outputs,
      t = sim.t,
    )
    calsim = calibration.CalibrationSim(
      spec['title'],
      sim = sim,
      targets = targets,
      verbose = True,
    )
    calsim.sim.params['tau'].update(0.05)
    opt = calsim.optimize(eps=0.1)
    calsim.fitted_params().save(fname_fit(calsim.name))

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
      'tpaf-WH': specfun(
        outputs = ['tpaf-WH'],
        selectors = ['all'],
        ylabel = 'tPAF of High Risk Women',
        ylim = [0,1]),
      'fit': specfun(
        fun = fitfun,
        outputs = ['prevalence'],
        selectors = [],
      ),
    # default: return all except 'fit', else, only user specified
    }.items() if (name in todo) or ((len(todo)==0) and not (name in ['fit','tpaf-fsw'] )) }

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
      title = variants.make_title(*variants.parse_name(varname)),
      save = fname_fig(plot,varname)
    )
    plt_close()

def run_sims(todo=None):
  todo = todo if todo is not None else []
  for i,(name,sim) in enumerate(variants.get_sims().items()):
    run_sim(sim,todo,name)
  
def print_turnover(name,sim):
  print('-'*50)
  print(name)
  print(sim.model.params['zeta'].islice(t=50,ki='M'))
  print(sim.model.params['dur'].islice(t=50,ki='M'))
  print(sim.model.params['pe'].islice(t=50,ki='M'))

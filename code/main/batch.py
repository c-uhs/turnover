import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
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

def fname_fig(plot,varname,ftype='pdf'):
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
    spec['fun'] = spec.pop('fun',plotfun)
    return spec

  def plotfun(**spec):
    name = spec.pop('name')
    sim.plot(**spec)

  def fitfun(**spec):
    # special case: fit model to targets, save fitted parameters (from params.json)
    # initialize the targets and calsim
    name = spec.pop('name')
    nu,G,Z = variants.parse_name(spec['name'])
    if G > 1:
      targets = system.get_targets(
        sim.model.spaces,
        sim.model.select,
        sim.outputs,
        t = sim.t,
      )
      calsim = calibration.CalibrationSim(
        spec['name'],
        sim = sim,
        targets = targets,
        verbose = True,
      )
      calsim.sim.params['tau'].update(0.05)
      opt = calsim.optimize(plot='tmp.pdf',ftol=1e-3)
      utils.savejson(fname_fit(calsim.name),calsim.fitted_params().todict())

  # dictionary of specifications for the plots
  specs = { name:spec for name,spec in \
    { 'N': specfun(
        output = 'N',
        selectors = ['high','med','low'],
        ylim = [0,5000]),
      'X-sir': specfun(
        output = 'X',
        selectors = ['S','I','R'],
        ylim = [0,1]),
      'X-groups': specfun(
        output = 'X',
        selectors = ['high','med','low'],
        ylim = [0,1]),
      'prevalence': specfun(
        output = 'prevalence',
        selectors = ['high','med','low'],
        ylim = [0,1]),
      'incidence': specfun(
        output = 'incidence',
        selectors = ['high','med','low'],
        ylim = [0,100]),
      'incidence-abs': specfun(
        output = 'incidence-abs',
        selectors = ['high','med','low'],
        ylim = [0,50]),
      'cum-infect': specfun(
        output = 'cum-infect',
        selectors = ['high','med','low'],
        ylim = [0,1000]),
      'tpaf-high': specfun(
        output = 'tpaf-high',
        selectors = ['all'],
        # ylabel = 'tPAF of High Risk Women',
        ylim = [0,1]),
      'fit': specfun(
        fun = fitfun,
        output = 'prevalence',
        selectors = []),
      # 'prop-from-WH': specfun(
      #   output = 'prop-from-WH',
      #   selectors = ['infected'],
      #   ylabel = 'testing ...',
      #   ylim = [0,1]),
    # default: return all except 'fit', else, only user specified
    }.items() if (name in todo) or ((len(todo)==0) and not (name in ['fit','tpaf-high'] )) }

  # initialize the required outputs
  outputs = utils.unique(utils.flatten(spec['output'] for spec in specs.values()))
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
      show  = False,
      name  = varname,
      title = variants.make_title(*variants.parse_name(varname)),
      save  = fname_fig(plot,varname)
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

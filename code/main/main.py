import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot()

import numpy as np
config.numpy()
from copy import deepcopy
from collections import OrderedDict as odict
import matplotlib.pyplot as plt

import system
import modelutils
import calibration
import utils

def get_sims_structure(outputs=[]):
  # build a dictionary of model variants, starting from the most complicated
  specs = system.get_specs()
  sims = odict()
  for nu in [specs['params']['mu'],specs['params']['nu']]:
    for G in [3,1]:
      for Z in [1,0]:
        # get the default model
        model = system.get_model()
        # define the growth rate
        model.params['nu'].update(nu)
        # define the number of of groups
        if G == 1:
          if Z != 0: continue # Z is irrelevant for G = 1, so only need 1 variant
          model.collapse(['ii'])
        # define turnover via specified durations
        if Z == 0:
          model.params['dur'].update(np.nan)
          model.params['zeta'].update(np.nan)
        # add simulation to dictionary
        sims.update({
          'nu={}_G={}_Z={}'.format(nu,G,Z):
          system.get_simulation(model,outputs=outputs)
        })
  return sims

def get_sims_zeta(outputs=[]):
  def get_zeta_specs(nu,mu,pe,dur,zeta):
    # build set of sufficient constraints given zeta with some entries zero:
    #   zeta will be either zero or nan (calculated)
    #   dur will be either input value of nan (calculated)
    #   pe will be nan (calculated) always
    # NOTE did not check validity of this constraint building for G > 3
    zeta = deepcopy(zeta)
    for i,zi in enumerate(zeta):
      zc = sum(np.isnan(zi))
      if zc == 0:
        # calculate dur if all zeta for this i specified
        dur[i] = np.nan
      if zc == 2:
        # specify zeta for this i if all to be calculated (even split)
        zsum = (1/dur[i] - mu)
        zeta[i] = [zsum/2 if (i != j) else (0) for j in range(3)]
      # calculate pe always
      pe[i] = np.nan
    return {'zeta':zeta,'dur':dur,'pe':pe}

  def zeta_str(zeta):
    return str(
      [
        ['.' if i==j else
         '0' if zij == 0 else
         'x' if np.isnan(zij) else
         'Z'
        for j,zij in enumerate(zi)]
      for i,zi in enumerate(zeta)]
    ).replace('\'','').replace(', ','')

  specs = system.get_specs()
  sims = odict()
  zeta = np.zeros((3,3))
  for i,ij in enumerate([None,(0,2),(0,1),(1,2),(2,0),(2,1),(1,0)]):
    model = system.get_model()
    if ij:
      zeta[ij[0],ij[1]] = np.nan
    for ki in ['M','W']:
      zspecs = get_zeta_specs(nu = model.params['nu'],
                              mu = model.params['mu'],
                              pe = model.params['pe'].islice(ki=ki),
                              dur = model.params['dur'].islice(ki=ki),
                              zeta = zeta)
      for param in ['zeta','dur','pe']:
        model.params[param].update(zspecs[param],ki=ki)
    sims.update({
      '({})-{}'.format(i,zeta_str(zspecs['zeta'])):
      system.get_simulation(model,outputs=outputs)
    })
  return sims

def runsim(sim,plots,variant,vset):
  # run a single simulation for a single model variant, and generate named plots
  specs = system.get_specs()

  def savename(vset,plot,variant,ftype='png'):
    # filename for the plot
    fname = plot+'_'+variant+'.'+ftype
    return os.path.join(config.path['root'],'outputs','figs','plots',vset,fname)

  def specfun(**spec):
    # clean-up the plot specifications
    spec['selectors'] = [sim.model.select[name] for name in spec['selectors']]
    spec['fun'] = spec.pop('fun',lambda **kwargs: sim.plot(**kwargs))
    return spec

  def tpafplotfun(**spec):
    # special plotting function for tPAF since it is not a model output
    spec.pop('outputs')
    pop = spec.pop('pop')
    tpaf = modelutils.tpaf(sim,pop,dxinf='xlam',beta='beta',copy=True)
    modelutils.plot(sim.t,tpaf,**spec)

  def fitfun(**spec):
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
      'tpaf-h': specfun(
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
    }.items() if (name in plots) or (len(plots)==0) }

  # initialize the required outputs
  outputs = utils.unique(utils.flatten(spec['outputs'] for spec in specs.values()))
  sim.init_outputs(system.get_outputs(
    spaces = sim.model.spaces,
    select = sim.model.select,
    t = sim.t,
    names = outputs))
  # solve the system
  print('> {}'.format(variant),flush=True)
  sim.solve()
  # generate and save the specified plots
  for plot,spec in specs.items():
    print('  + {}'.format(plot),flush=True)
    spec.pop('fun')(**spec,
      show = False,
      title = variant,
      save = savename(vset,plot,variant))
    plt.close()

if __name__ == '__main__':

  # variants w.r.t. model structure
  # for variant,sim in get_sims_structure(outputs=[]).items():
  #   runsim(sim,[],variant,'structure')

  # variants w.r.t. values of zeta
  for i,(variant,sim) in enumerate(get_sims_zeta(outputs=[]).items()):
    runsim(sim,['fit'],variant,'zeta')

    # # double check key turnover parameters after solving
    # print('-'*50)
    # print(variant)
    # print(sim.model.params['zeta'].islice(t=2000,ki='M'))
    # print(sim.model.params['dur'].islice(t=2000,ki='M'))
    # print(sim.model.params['pe'].islice(t=2000,ki='M'))
      # break


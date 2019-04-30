import os,sys;
sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config'));
import config
config.epimodel()
# epi-model imports
from utils import atleast,xdict,xfilter
from space import Dimension,Space,Array
from simulation import Model,Simulation,ParameterSet
from elements import Selector
import modelutils
import outpututils
import initutils
import transmit
import mixing
# external imports
import numpy as np

def dxfun(X,t,P,dxout={}):
  # TODO: validate this!
  def foifun(Xi,Ci,beta):
    Xp  = modelutils.partner(X)
    Cp  = modelutils.partner(Ci)
    XC  = Xp.isum(['hp']) * Cp
    rho = XC / XC.isum(['ip'])
    rho[np.isnan(rho)] = 1
    rho = rho.expand(beta.space.subspace(['hp'],keep=False))
    return Ci * (
      rho * (beta * Xp).isum(['hp']) / Xp.isum(['hp'])
    ).isum(['ip'])
    # / TODO
  # dxfun
  dXi = X*0 # entry to  all compartments
  dXo = X*0 # exit from all compartments
  # births and deaths
  dXi.update(P['nu']*P['pe']*X.sum(), hi=['S'], accum=np.add)
  dXo.update(P['mu']*X, accum=np.add)
  # turnover: ii -> ip ("from" group is index, not "to")
  if X.space.dim('ii').n > 1:
    Xi = X.expand(Space(X.space.dims+[modelutils.partner(X.space.dim('ii'))]),norm=False)
    XZ = Xi * atleast(P['zeta'],3,end=+1,cp=True)
    dXi.update(XZ.isum(['ii']),accum=np.add)
    dXo.update(XZ.isum(['ip']),accum=np.add)
  # force of infection: S -> I
  lam  = foifun(X,P['C'],P['beta'])
  xlam = lam.iselect(hi=['S'])*X.iselect(hi=['S'])
  dXi.update(xlam, hi=['I'], accum=np.add)
  dXo.update(xlam, hi=['S'], accum=np.add)
  # treatment: I -> T
  treat = X.iselect(hi=['I']) * P['tau']
  dXi.update(treat, hi=['R'], accum=np.add)
  dXo.update(treat, hi=['I'], accum=np.add)
  # dxout
  if dxout: # TODO: this should be a decorator
    lvars = locals()
    for v in dxout:
      if v in lvars:
        dxout.update({v:locals()[v]})
  # return
  return dXi - dXo

def initfun(model):
  P = model.params
  # define beta
  P['beta'].update(0)
  P['beta'].update(P['beta-i'],hp=['I'])
  # turnover
  if P['zeta'].space.dim('ii').n > 1:
    if np.all(np.isnan(P['dur'])): # TODO: this is not a robust condition
      P['zeta'].update(0)
    else:
      turnover = transmit.turnover(
        nu   = P['nu'],
        mu   = P['mu'],
        px   = P['px'],
        pe   = P['pe'],
        zeta = P['zeta'],
        dur  = P['dur'],
        warn = True,
      )
      P['zeta'].update(turnover['zeta'])
      P['dur'].update(turnover['dur'])
      P['pe'].update(turnover['pe'])
  # initial condition
  model.X0.update(P['N0']*P['px'],hi=['S'])

def infectfun(sim,N=1):
  sim.init_x(transmit.transfer(
    X    = sim.X,
    src  = {'hi':['S']},
    dst  = {'hi':['I']},
    both = {'t':[sim.t[0]]},
    N    = N,
  ))
  return sim

def get_targets(spaces,select,outputs,t=None,names=None):
  specdir = os.path.join(config.path['root'],'code','main','specs')
  targets = initutils.objs_from_json(
    initutils.make_target,
    os.path.join(specdir,'targets.json'),
    space   = spaces['super'],
    select  = select,
    outputs = outputs,
  )
  names = [target.name for target in targets.values()] if names is None else flatten(names)
  return xdict(xfilter(targets.values(),name=names),ordered=True)

def get_outputs(spaces,select,t,names=None,**kwargs):
  specs = {
    'N':             {},
    'X':             {},
    'prevalence':    {'si':'infected'},
    'incidence':     {'ss':'S'},
    'incidence-abs': {'ss':'S'},
    'cum-infect':    {'ss':'S'},
    'tpaf-high':     {'beta':'beta','ss':'S','mode':'to'},
    'inf-ratio':     {'ss':'S','si':'infected'},
  }
  return xdict([
    outpututils.make_output(name,t=t,spaces=spaces,**specs[name])
    for name in names
  ])

def get_specs():
  specdir = os.path.join(config.path['root'],'code','main','specs')
  dims   = initutils.objs_from_json(Dimension,            os.path.join(specdir,'dimensions.json'))
  spaces = initutils.objs_from_json(initutils.make_space, os.path.join(specdir,'spaces.json'),dims=dims.values())
  params = initutils.objs_from_json(initutils.make_param, os.path.join(specdir,'params.json'),space=spaces['super'])
  select = initutils.objs_from_json(Selector,             os.path.join(specdir,'selectors.json'))
  return {
    'dims':   dims,
    'spaces': spaces,
    'params': params,
    'select': select,
  }

def get_model():
  specs = get_specs()
  return Model(
     X0      = Array(0,specs['spaces']['index']),
     dxfun   = dxfun,
     spaces  = specs['spaces'],
     select  = specs['select'],
     params  = ParameterSet(specs['params']),
     initfun = initfun,
   )

def get_t(dt=0.5,tmin=0,tmax=200):
  return np.around(np.arange(tmin, tmax+1e-6, dt), 6)

def get_simulation(model,infect=True,outputs=[],t=None):
  infect = atleast(model.params['infect'],2) if 'infect' in model.params else \
           int(infect)
  t = t if t is not None else get_t()
  outputs = get_outputs(model.spaces,model.select,t=t,names=outputs)
  sim = Simulation(model,t,outputs=outputs)
  return infectfun(sim,N=infect)

import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
# epi-model imports
from space import *
from simulation import *
import outputfuns
import initutils
import transmit
import test
# external imports
import numpy as np

def dxfun(X,t,P,dxout={}):
  dX = X*0
  # births and deaths
  dX.update(P['nu']*P['pe']*X.sum(), hi='S', accum=np.add)
  dX.update(P['mu']*X, accum=np.subtract)
  # turnover
  if P['zeta'].space.dim('ii').n > 1:
    Xi = X.expand(Space(X.space.dims+[modelutils.partner(X.space.dim('ii'))]),norm=False)
    for ki in ['M','W']:
      XZk = Xi.iselect(ki=ki) * P['zeta'].iselect(ki=ki)
      dX.update(XZk.isum('ip'),ki=ki,accum=np.subtract)
      dX.update(XZk.isum('ii'),ki=ki,accum=np.add)
  # force of infection
  lamp = transmit.lambda_fun(X,P['C'],P['beta'].islice(t=t),P['eps'],dxout)
  xlam = lamp.isum('p').iselect(hi='S')*X.iselect(hi='S')
  transmit.transfer(dX, src={'hi':'S'}, dst={'hi':'I'}, N = xlam)
  # treatment
  transmit.transfer(dX, src={'hi':'I'}, dst={'hi':'T'}, N = X.iselect(hi='I') * P['tau'])
  # # full recovery -> for SIT model, we have no recovery
  # transmit.transfer(dX, src={'hi':'T'}, dst={'hi':'S'}, N = X.iselect(hi='T') * P['gamma'])
  # dxout
  dxout.update({v:locals()[v] for v in dxout if v in locals()})
  # return
  return dX

def initfun(model,spaces):
  P = model.params
  # beta
  sbeta = spaces['super'].union(P['f-beta-h'].space)
  P['beta'] = P['beta'].expand(sbeta) * P['f-beta-h'].expand(sbeta)
  # turnover
  if P['zeta'].space.dim('ii').n > 1:
    for ki in ['M','W']:
      turnover = transmit.turnover(nu = P['nu'],
                                   mu = P['mu'],
                                   px = P['px'].islice(ki=ki),
                                   pe = P['pe'].islice(ki=ki),
                                   zeta = P['zeta'].islice(ki=ki),
                                   dur = P['dur'].islice(ki=ki),
                                   warn = True)
      P['zeta'].update(turnover['zeta'],ki=ki)
      P['dur'].update(turnover['dur'],ki=ki)
      P['pe'].update(turnover['pe'],ki=ki)
  # initial condition
  model.X0.update(P['N0']*P['px'],hi='S')

def infectfun(sim,N=1):
  sim.init_x(transmit.transfer(
    X = sim.X,
    src = {'hi':'S'},
    dst = {'hi':'I'},
    both = {'t':sim.t[0]},
    N = N))
  return sim

def get_outputs(spaces,select,t=None,names=None):
  def tspace(space):
    return space.union(Space([modelutils.tdim(flatten(t))]))
  outputs = [
    Output('N',
           space = tspace(spaces['index']),
           fun = lambda sim: sim.X,
           accum = np.sum,
           wax = False),
    Output('X',
           space = tspace(spaces['index']),
           fun = lambda sim: sim.X / sim.X.isum('t',keep=True),
           accum = np.sum,
           wax = False),
    Output('prevalence',
           space = tspace(spaces['index']),
           fun = lambda sim: outputfuns.prevalence(sim),
           accum = np.average,
           wax = True),
    Output('incidence',
          space = tspace(spaces['super'].subspace(['ki','ii','p'],keep=True)),
          fun = lambda sim,t,lamp: outputfuns.incidence(sim,lam=lamp,t=t,per=1000,ss='S'),
          accum = np.average,
          wax = True,
          calc = 'peri', dxout = ['lamp']),
    Output('incidence-abs',
          space = tspace(spaces['super'].subspace(['ki','ii','p'],keep=True)),
          fun = lambda sim,t,lamp: outputfuns.incidence(sim,lam=lamp,t=t,per=False,ss='S'),
          accum = np.average,
          wax = True,
          calc = 'peri', dxout = ['lamp']),
    Output('cum-infect',
          space = tspace(spaces['super'].subspace(['ki','ii','p'],keep=True)),
          fun = lambda sim,t,lamp: outputfuns.incidence(sim,lam=lamp,t=t,per=False,ss='S'),
          accum = np.sum,
          wax = False,
          calc = 'peri', dxout = ['lamp'], cum = True),
  ]
  names = [output.name for output in outputs] if names is None else flatten(names)
  return xdict(xfilter(outputs,name=names))

def get_specs():
  specdir = os.path.join(config.path['root'],'code','main','specs')
  dims   = initutils.objs_from_json(Dimension,            os.path.join(specdir,'dimensions.json'))
  spaces = initutils.objs_from_json(initutils.make_space, os.path.join(specdir,'spaces.json'),dims=dims.values())
  params = initutils.objs_from_json(initutils.make_param, os.path.join(specdir,'params.json'),space=spaces['super'])
  accum  = initutils.objs_from_json(initutils.make_accum, os.path.join(specdir,'accumulators.json'),params=params)
  select = initutils.objs_from_json(Selector,             os.path.join(specdir,'selectors.json'))
  return {
    'dims': dims,
    'spaces': spaces,
    'params': params,
    'select': select,
    'accum': accum,
  }

def get_model():
  specs = get_specs()
  return Model(X0 = Array(0,specs['spaces']['index']),
               dxfun = dxfun,
               spaces = specs['spaces'],
               select = specs['select'],
               params = ParameterSet(specs['params'],accum=specs['accum']),
               initfun = lambda model: initfun(model,model.spaces))

def get_simulation(model,infect=True,outputs=[]):
  infect = atleast(model.params['infect'],3) if 'infect' in model.params else \
           int(infect)
  t = np.around(np.arange(1975, 2025+1e-6, 0.5),6)
  outputs = get_outputs(model.spaces,model.select,t=t,names=outputs)
  sim = Simulation(model,t,outputs=outputs)
  return infectfun(sim,N=infect)

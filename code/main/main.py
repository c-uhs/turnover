import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot()
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
  Xi = X.expand(Space(X.space.dims+[modelutils.partner(X.space.dim('ii'))]),norm=False)
  for ki in ['M','W']:
    XZk = Xi.iselect(ki=ki) * P['zeta'].iselect(ki=ki)
    dX.update(XZk.isum('ip'),ki=ki,accum=np.subtract)
    dX.update(XZk.isum('ii'),ki=ki,accum=np.add)
  # force of infection
  lamp = transmit.lambda_fun(X,P['C'],P['beta'].islice(t=t),P['eps'],dxout)
  transmit.transfer(dX, src={'hi':'S'}, dst={'hi':'I'}, N = lamp.isum('p').iselect(hi='S')*X.iselect(hi='S'))
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
  for ki in ['M','W']:
    P['zeta'].update(transmit.zeta_fun( P['nu'], P['mu'],
        P['pe'].islice(ki=ki), P['px'].islice(ki=ki), P['dur'].islice(ki=ki)),
      ki=ki)
  # initial condition
  model.X0.update(P['N0']*P['px'],hi='S')

def infectfun(sim):
  sim.init_x(transmit.transfer(
    X = sim.X,
    src = {'hi':'S'},
    dst = {'hi':'I'},
    both = {'t':sim.t[0]},
    N = 1))
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
          fun = lambda sim,t,lamp: outputfuns.incidence(sim,lam=lamp,t=t,per=1000,s='S'),
          accum = np.average,
          wax = True,
          calc = 'peri', dxout = ['lamp']),
    Output('cum-infect',
          space = tspace(spaces['super'].subspace(['ki','ii','p'],keep=True)),
          fun = lambda sim,t,lamp: outputfuns.incidence(sim,lam=lamp,t=t,per=False),
          accum = np.sum,
          wax = False,
          calc = 'peri', dxout = ['lamp'], cum = True),
  ]
  names = [output.name for output in outputs] if names is None else flatten(names)
  return xdict(xfilter(outputs,name=names))

specdir = os.path.join(config.path['root'],'code','main','specs')
dims   = initutils.objs_from_json(Dimension,            os.path.join(specdir,'dimensions.json'))
spaces = initutils.objs_from_json(initutils.make_space, os.path.join(specdir,'spaces.json'),dims=dims.values())
params = initutils.objs_from_json(initutils.make_param, os.path.join(specdir,'params.json'),space=spaces['super'])
select = initutils.objs_from_json(Selector,             os.path.join(specdir,'selectors.json'))
accum  = initutils.objs_from_json(initutils.make_accum, os.path.join(specdir,'accumulators.json'),params=params)
model = Model(X0 = Array(0,spaces['index']),
             dxfun = dxfun,
             params = params,
             select = select,
             initfun = lambda model: initfun(model,spaces))

# model.params.collapse(accum,idxs=['ii'])
# model.params['zeta']

t = np.around(np.arange(1975, 2025+1e-6, 0.1),6)
outputs = get_outputs(spaces,select,t=t,names=['N','X','prevalence'])
sim = infectfun(Simulation(model,t,outputs=outputs))

# # 1-off plotting
# sim.solve()
# sim.plot(outputs=['X'],selectors=[model.select[name] for name in ['S','I','T']])
# sim.plot(outputs=['incidence'],selectors=[model.select[name] for name in ['WH','MH','WM','MM','WL','ML']])

# # gif plotting
def gifplot(dur):
  # sim.model.params['nu'].update(nu)
  sim.model.params['dur'].update(dur,ki='W',ii='H')
  sim.init_model(sim.model)
  sim.init_params()
  sim.solve()
  sim.plot(outputs=['prevalence'],
           # selectors=[model.select[name] for name in ['S','I','T']],
           selectors=[model.select[name] for name in ['WH','MH','WM','MM','WL','ML']],
           show=False)
  import matplotlib.pyplot as plt
  plt.gca().set_ylim([0,1])

# using decorator from utils
gif(fname='test.gif',clf=True,verbose=True,overlay=(0.05,0.95),
    # nu=np.linspace(0.03,0.06,5),
    dur=np.linspace(5,15,5),
    # eps=np.linspace(0,1,6),
  )(gifplot)
    
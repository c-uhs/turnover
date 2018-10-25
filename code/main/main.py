import os,sys; sys.path.append(os.path.join((lambda r,f:f[0:f.index(r)+len(r)])('code',os.path.abspath(__file__)),'config')); import config
config.epimodel()
config.plot()
# epi-model imports
from space import *
from simulation import *
import initutils
import transmit
# external imports
import numpy as np

ATTRDEATHS = False

def dxfun(X,t,P,dxout={}):
  dX = X*0
  # births and deaths
  dX.update(P['nu']*P['pe']*X.sum(), hi='S', accum=np.add)
  dX.update(P['mu']*X, accum=np.subtract)
  # attributable death
  if ATTRDEATHS:
    dX.update(X.iselect(hi='I')*P['phi'], hi='I', accum=np.subtract)
  # turnover
  Xi = X.expand(Space(X.space.dims+[modelutils.partner(X.space.dim('ii'))]),norm=False)
  for ki in ['M','W']:
    phi = P['phi'] * X.islice(ki=ki,hi='I') / X.islice(ki=ki).isum('hi')
    if not ATTRDEATHS:
      phi = 0*phi
    zeta = transmit.zeta_fun(P['nu'],P['mu'],phi,P['pe'].islice(ki=ki),P['px'].islice(ki=ki))
    XZk = Xi.iselect(ki=ki) * zeta
    dX.update(XZk.isum('ip'),ki=ki,accum=np.subtract)
    dX.update(XZk.isum('ii'),ki=ki,accum=np.add)
  # force of infection
  lam = transmit.lambda_fun(X,P['C'],P['beta'],P['eps'],dxout).isum('p')
  transmit.transfer(dX, src={'hi':'S'}, dst={'hi':'I'}, N = lam.iselect(hi='S')*X.iselect(hi='S'))
  # treatment
  transmit.transfer(dX, src={'hi':'I'}, dst={'hi':'R'}, N = X.iselect(hi='I') * P['tau'])
  # loss of immunity
  transmit.transfer(dX, src={'hi':'R'}, dst={'hi':'S'}, N = X.iselect(hi='R') * P['gamma'])
  # dxout
  dxout.update({v:locals()[v] for v in dxout if v in locals()})
  # return
  return dX

def initfun(model,spaces):
  P = model.params
  P['beta'] = P['beta'].expand(spaces['super']) * P['f-beta-h'].expand(spaces['super'])
  model.X0.update(P['N0']*P['px'],hi='S') # initial condition

def infectfun(sim):
  sim.init_x(transmit.transfer(
    X = sim.X,
    src = {'hi':'S'},
    dst = {'hi':'I'},
    both = {'t':sim.t[0]},
    N = 1))
  return sim

specdir = os.path.join(config.path['root'],'code','main','specs')
dims   = initutils.objs_from_json(Dimension,            os.path.join(specdir,'dimensions.json'))
spaces = initutils.objs_from_json(initutils.make_space, os.path.join(specdir,'spaces.json'),dims=dims.values())
select = initutils.objs_from_json(Selector,             os.path.join(specdir,'selectors.json'))
params = initutils.objs_from_json(initutils.make_param, os.path.join(specdir,'params.json'),space=spaces['super'])
model = Model(X0 = Array(0,spaces['index']),
             dxfun = dxfun,
             params = params,
             initfun = lambda model: initfun(model,spaces))
# sim = model.equilibriate()
t = np.around(np.arange(1975, 2025+1e-6, 1),6)
sim = infectfun(Simulation(model,t))
sim.solve()
sim.X *= (sim.X.isum('t',keep=True))**-1
# sim.plot(selectors=[select[name] for name in ['S','I','R']])
sim.plot(selectors=[select[name] for name in ['WH','MH','WM','MM','WL','ML']])

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

def dxfun(X,t,P,dxout={}):
  dX = X*0
  # force of infection
  lam = transmit.lambda_fun(X,P['C'],P['beta'],P['eps'],dxout).isum('p')
  transmit.transfer(dX, src={'hi':'S'}, dst={'hi':'I'},
    N = lam.iselect(hi='S')*X.iselect(hi='S'))
  # treatment
  transmit.transfer(dX, src={'hi':'I'}, dst={'hi':'R'},
    N = X.iselect(hi='I') * P['tau'])
  # loss of immunity
  transmit.transfer(dX, src={'hi':'R'}, dst={'hi':'S'},
    N = X.iselect(hi='R') * P['gamma'])
  # locals
  dxout.update({v:locals()[v] for v in dxout if v in locals()})
  # return
  return dX

def initfun(model,spaces):
  P = model.params
  model.X0.update(P['N0']*P['px'],hi='S')
  P['beta'] = P['beta'].expand(spaces['super'])

def infectfun(sim):
  sim.init_x(transmit.transfer(
    X = sim.X.islice(t=sim.t[0]),
    src = {'hi':'S'},
    dst = {'hi':'I'},
    N   = 1))
  return sim

dims   = initutils.objs_from_json(Dimension,            os.path.join('specs','dimensions.json'))
spaces = initutils.objs_from_json(initutils.make_space, os.path.join('specs','spaces.json'),dims=dims.values())
select = initutils.objs_from_json(Selector,             os.path.join('specs','selectors.json'))
params = initutils.objs_from_json(initutils.make_param, os.path.join('specs','params.json'),space=spaces['super'])
model = Model(X0 = Array(0,spaces['index']),
             dxfun = dxfun,
             params = params,
             initfun = lambda model: initfun(model,spaces))
# sim = model.equilibriate()
t = np.around(np.arange(1975, 2025+1e-6, 0.1),6)
sim = infectfun(Simulation(model,t))
sim.solve()
sim.plot(selectors=[select[name] for name in ['S','I','R']])
sim.plot(selectors=[select[name] for name in ['WH','MH','WM','MM','WL','ML']])

# P = np.array([0.4,0.3,0.2,0.1])
# print(zetafun(0.3,0.2,P,P,[4,np.nan,np.nan,np.nan]))
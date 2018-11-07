from scipy.optimize import minimize
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np

xi  = np.array([0.10,0.20,0.70])
xo  = np.array([0.10,0.20,0.70])
ri  = 0.05
ro  = 0.02
dur = np.array([5,np.nan,np.nan])
Aprime = None # [1,0,0,0,0,0]
bprime = None # [0.05]
L2 = 0.00

G = len(xo)
# define the LHS
b = ri*(xo - xi).T
# vectors of indices for convenience
iz = [(i,j) for i in range(G) for j in range(G) if i is not j]
iv = [i*G+j for i in range(G) for j in range(G) if i is not j]
# define the RHS system matrix
A = np.array(\
  [ [-xo[zi[0]] if zi[0] == i else 
      xo[zi[0]] if zi[1] == i else 0
      for zi in iz
    ] for i in range(G)
  ])
# append duration-based constraints
if dur is not None:
  for id,di in enumerate(dur):
    if np.isfinite(di):
      b = np.concatenate((b, [(di**(-1)) - ro]))
      A = np.concatenate((A, [[1 if (i == id) else 0
                              for i in range(G)
                              for j in range(G) if i is not j]]),axis=0)
# append any additional constraints
if (bprime is not None) and (Aprime is not None):
  b = np.concatenate(( b, np.atleast_1d(bprime) ), axis=0)
  A = np.concatenate(( A, np.atleast_2d(Aprime) ), axis=0)
# solve the system
eps = 1e-8
jfun = lambda z: np.linalg.norm((np.dot(A,z)-b),2) + L2*np.linalg.norm(z,2)
z0 = np.zeros((G*G-G,1))
out  = minimize(jfun, z0, bounds = [(0.00,0.50) for i in z0],
                          method = 'TNC', options = {'ftol':eps})
print(out['fun'])
z = out['x']
# return zeta as a matrix
zeta = np.array(\
  [ [ z[iv.index(i*G+j)] if i*G+j in iv else 0
      for j in range(G)
    ] for i in range(G)
  ])
print(zeta)

N = len(z)
for ij in combinations(range(N),2):
  zmax = max(max(z[ij[0]],z[ij[1]]),1e-6)
  zi,dzi = np.linspace(0,2*zmax,20,retstep=True)
  zj,dzj = np.linspace(0,2*zmax,20,retstep=True)
  J = np.array([[ jfun([ i if k == ij[0] else  j if k == ij[1] else zk \
    for k,zk in enumerate(z)]) \
      for i in zi]
        for j in zj])

  plt.subplot(N-1,N-1,ij[0]*(N-1)+ij[1])
  plt.imshow(-np.log(J+2*eps),cmap='magma',interpolation='bilinear')
  zix,zjx = z[ij[0]], z[ij[1]]
  plt.plot(zix/dzi, zjx/dzj, 'bx', markersize=10)
  title = 'z{}{} vs z{}{}'.format(\
    *[e+1 for e in iz[ij[0]]],
    *[e+1 for e in iz[ij[1]]])
  plt.title(title)
  print(title+': ({:.6f},{:.6f})'.format(zix,zjx),flush=True)
plt.show()

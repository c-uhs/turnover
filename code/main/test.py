from scipy.optimize import minimize
import numpy as np

WARN = False

def zeta_fun(ri,ro,xi,xo,D=None,bprime=None,Aprime=None):
  r"""Compute the turnover matrix to yield a steady-state population distribution

  Solve for the internal transitions matrix $\\zeta$ required
  to maintain constant population proportions $\\bm{{x}}$
  for given entry and exit rates and distributions.
  Depending on the length of $\\bm{{x}}$, the system may be underdetermined;
  if so, ${{\\left\|\\left\|\\zeta\\right\|\\right\|}}\_2$ is also minimized.
  Additional constraints of the form $b' = A'\\bm{{z}}$
  can also be provided to help constrain the solution.

  .. include:: static/eq/zeta.rst

  Args:
    ri: $\\nu$ entry rate
    ro: $\\mu$ exit rate
    xi: $\\bm{{x}}_i$ entry population distrubution (vector)
    xo: $\\bm{{x}}$ target population distribution (vector)
    D: $\\bm{{D}}$ target duration in each group (vector)
    bprime: $b'$ additional constraints LHS
    Aprime: $A'$ additional rows of A (RHS)

  Returns:
    the transitions matrix $\\zeta$
  """
  N = len(xo)
  # define the LHS
  b = ri*(xo - xi).T
  # vectors of indices for convenience
  iz = [(i,j) for i in range(N) for j in range(N) if i is not j]
  iv = [i*N+j for i in range(N) for j in range(N) if i is not j]
  # define the RHS system matrix
  A = np.array(\
    [ [-xo[zi[0]] if zi[0] == i else 
        xo[zi[0]] if zi[1] == i else 0
        for zi in iz
      ] for i in range(N)
    ])
  # append duration-based constraints
  if D is not None:
    for id,di in enumerate(D):
      if np.isfinite(di):
        b = np.concatenate((b, [(di**(-1)) - ro]))
        A = np.concatenate((A, [[1 if (i == id) else 0
                                for i in range(N)
                                for j in range(N) if i is not j]]),axis=0)
  # append any additional constraints
  if (bprime is not None) and (Aprime is not None):
    b = np.concatenate(( b, np.atleast_1d(bprime) ), axis=0)
    A = np.concatenate(( A, np.atleast_2d(Aprime) ), axis=0)
  # solve the system
  jfun = lambda z: np.linalg.norm((np.dot(A,z)-b),2)
  # z0   = np.dot(np.linalg.pinv(A),b)
  z0 = np.zeros((N*N-N,1))
  out  = minimize(jfun, z0, bounds = [(0.00,0.50) for i in z0],
                            method = 'TNC', options = {'ftol':1e-8})
  z  = out['x']
  if WARN:
    if not (out['fun'] < 1e-8):
      print('Warning: zeta_fun did not reach zero; final jfun = {}'.format(out['fun']))
    if not out['success']:
      print('Warning: zeta_fun did not converge; final jfun = {}'.format(out['fun']))
  # return zeta as a matrix
  zeta = np.array(\
    [ [ z[iv.index(i*N+j)] if i*N+j in iv else 0
        for j in range(N)
      ] for i in range(N)
    ])
  return zeta
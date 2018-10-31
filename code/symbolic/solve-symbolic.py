# README
# Goal: Solve for the values of the turnover matrix (z11), z12, z13, ...
# in terms of arbitrary distribution x1, x2, x3, ... and rates of change b1, b2, b3, ...
# Status: seems possible but sympy gives no solution
# Author: Jesse Knight, Oct 2018

import sympy as sp

N = 3
x = sp.Matrix([sp.Symbol('x'+str(i+1)) for i in range(N)])
z = [sp.Symbol('z'+str(i+1)+str(j+1)) for i in range(N) for j in range(N) if not i==j]
b = [sp.Symbol('b'+str(i+1)) for i in range(N)]

iz = [(i,j) for i in range(N) for j in range(N) if i is not j]
iv = [i*N+j for i in range(N) for j in range(N) if i is not j]

A = [ [-x[zi[0]] if zi[0] == i else 
        x[zi[0]] if zi[1] == i else 0
        for zi in iz
      ] for i in range(N)
    ]
    
Az  = sp.Matrix([[ai*zi for ai,zi in zip(row,z)] for row in A])
bAz = sp.Matrix([list(Az[ri,:])+[bi] for ri,bi in zip(range(N),b)])
print('SYSTEM:')
print('\n'.join(['0 = '+' + '.join([str(bAz[i,j]) for j in range(N*(N-1)+1)])  for i in range(N)]))
print('SOLUTION: [sympy.solve_linear_system]:')
print(sp.solve_linear_system(bAz,*z))

# Apparently no solution, but this doesn't seem right.
# The sympy package seems weak.
# Consider: https://octave.sourceforge.io/symbolic/

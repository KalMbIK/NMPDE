import numpy as np
from math import exp, sin, pi, fabs
import matplotlib.pyplot as plt

def solver(u0, M_, nu_, s_):
    u = np.array(u0)
    v = np.zeros_like(u)
    N = len(u0)
    nu_2 = nu_/2.
    s_2 = s_/2.
    sigma = [-nu_2+s_2,1.-s_,nu_2+s_2]
    for n in range(M_):
        v[0] = sigma[0]*u[1]+sigma[1]*u[0]+sigma[2]*u[N-1]
        v[1:-1] = sigma[0]*u[2:]+sigma[1]*u[1:-1]+sigma[2]*u[:-2]
        v[N-1] = sigma[0]*u[0]+sigma[1]*u[N-1]+sigma[2]*u[N-2]
        # v[0] = u[0]-nu_2*(u[1]-u[N-1])+s_2*(u[1]-2.*u[0]+u[N-1])
        # v[1:-1] = u[1:-1]-nu_2*(u[2:]-u[:-2])+s_2*(u[2:]-2.*u[1:-1]+u[:-2])
        # v[N-1] = u[N-1]-nu_2*(u[0]-u[N-2])+s_2*(u[0]-2.*u[N-1]+u[N-2])
        t = u
        u = v
        v = t
    return u

nu = 0.9
a = 1.

initialCondition1 = lambda x: sin(2.*pi*x)
initialCondition2 = lambda x: 1 if x < 0.5 else 0
exactSolution = lambda f, x, t: f(x-a*t)
maxNorm = lambda u, exact: np.fabs(u-exact).max()

upwind = lambda u0, M_: solver(u0, M_, nu, fabs(nu))
laxWendroff = lambda u0, M_: solver(u0, M_, nu, nu ** 2)

def problemSetup(solver, initCond, type, N, nu, a):
    X,h = np.linspace(0.,1.,N,endpoint=False,retstep=True)
    tau = nu*h/a
    M = int(N/nu+1)
    T = M*tau
    mu = fabs(a)*h/2.*(1-fabs(nu))

    exactSolution1a = lambda x,t: exactSolution(initCond,x,t)
    exactSolution1b = lambda x,t: exactSolution1a(x,t)*exp(-mu*4.*t*pi**2)
    sol = exactSolution1a if type == 'a' else exactSolution1b

    u0 = np.array([initCond(x) for x in X])
    exact = np.array([sol(x,T) for x in X])

    u = solver(u0, M)
    return maxNorm(u,exact)

for scheme, name in zip([upwind, laxWendroff],['upwind', 'LaxWendroff']):
    grids = [10**i for i in range(1, 6)]
    errors = []
    for n in grids:
        t = problemSetup(scheme,initialCondition1,'a',n,nu,a)
        errors.append(t)
    plt.loglog(grids, errors, label=name)

plt.legend(loc='best')
plt.show()
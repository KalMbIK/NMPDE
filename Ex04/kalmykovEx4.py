# import numpy as np
from math import exp, sin, pi, fabs, ceil, modf
import matplotlib.pyplot as plt

def my_zeros(n):
    return [0.]*n

def my_zeros_like(u):
    return my_zeros(len(u))

def my_array_copy(u):
    uu = my_zeros_like(u)
    for i in range(len(uu)):
        uu[i]=u[i]
    return uu

def my_vector_add(u1, u2):
    if len(u1) != len(u2):
        return None
    uu = my_zeros_like(u1)
    for i in range(len(uu)):
        uu[i] = u1[i]+u2[i]
    return uu

def my_vector_float_mult(alpha, u):
    uu = my_zeros_like(u)
    for i in range(len(uu)):
        uu[i] = alpha*u[i]
    return uu

def my_vector_sub(u1,u2):
    return my_vector_add(u1,my_vector_float_mult(-1., u2))

def my_vector_vectorize(f, u):
    uu = my_zeros_like(u)
    for i in range(len(u)):
        uu[i] = f(u[i])
    return uu

def my_vector_max_element(u):
    max = 0.
    for i in range(len(u)):
        if u[i]>= max:
            max = u[i]
    return max

def my_vector_linspace(init, end, npoints, endpoint ,retstep):
    if endpoint==True:
        step = (end-init)/(npoints-1)
    else:
        step = (end-init)/(npoints)
    grid = my_zeros(npoints)
    grid[0] = init
    for i in range(1, npoints):
        grid[i]=grid[i-1]+step
    if retstep==True:
        return grid, step
    else:
        return grid

def my_sign(x):
    if x < 0:
        return -1.
    elif x == 0:
        return 0.
    elif x > 0:
        return 1.

# Finite Difference Scheme
def solver(u0, M_, nu_, s_):
    u = my_array_copy(u0)
    v = my_zeros_like(u)
    # u = np.array(u0)
    # v = np.zeros_like(u)
    N = len(u0)
    nu_2 = nu_/2.
    s_2 = s_/2.
    sigma = [-nu_2+s_2,1.-s_,nu_2+s_2]
    for n in range(M_):
        v[0] = sigma[0]*u[1]+sigma[1]*u[0]+sigma[2]*u[N-1]
        for i in range(1, N-1):
            v[i] = sigma[0] * u[i+1] + sigma[1] * u[i] + sigma[2] * u[i-1]
            pass
        # v[1:-1] = sigma[0]*u[2:]+sigma[1]*u[1:-1]+sigma[2]*u[:-2]
        v[N-1] = sigma[0]*u[0]+sigma[1]*u[N-1]+sigma[2]*u[N-2]
        # v[0] = u[0]-nu_2*(u[1]-u[N-1])+s_2*(u[1]-2.*u[0]+u[N-1])
        # v[1:-1] = u[1:-1]-nu_2*(u[2:]-u[:-2])+s_2*(u[2:]-2.*u[1:-1]+u[:-2])
        # v[N-1] = u[N-1]-nu_2*(u[0]-u[N-2])+s_2*(u[0]-2.*u[N-1]+u[N-2])
        t = u
        u = v
        v = t
    return u

# Problem parameters
nu = 0.9
a = 1.

# Define initial conditions
initialCondition1 = lambda x: sin(2.*pi*x)
def initialCondition2(x):
    t = modf(x)
    if my_sign(x) >= 0:
    # if np.sign(x) >= 0:
        if t[0] >= 0.5:
            return 0
        else:
            return 1
    else:
        if fabs(t[0]) >= 0.5:
            return 1
        else:
            return 0

# Define some useful methods to evaluate exact solution, norm and schemes
exactSolution = lambda f, x, t: f(x-a*t)
maxNorm = lambda u, exact: my_vector_max_element(my_vector_vectorize(fabs,my_vector_sub(u,exact)))
# maxNorm = lambda u, exact: np.fabs(u-exact).max()

upwind = lambda u0, M_: solver(u0, M_, nu, fabs(nu))
laxWendroff = lambda u0, M_: solver(u0, M_, nu, nu ** 2)

# Input parametes:
# @solver - desirable scheme
# @initCond - initial conditions
# @type - 'a' - means that we are solving original advection equation,
# 'b' - means that we are solving modified equation for upwind scheme
# @N - number of nodes

def problemSetup(solver, initCond, type, N, nu, a):
    X,h = my_vector_linspace(0.,1.,N,endpoint=False,retstep=True)
    # X,h = np.linspace(0.,1.,N,endpoint=False,retstep=True)
    tau = nu*h/a
    M = int(ceil(N/nu))
    T = M*tau
    mu = fabs(a)*h/2.*(1-fabs(nu))

    exactSolution1a = lambda x,t: exactSolution(initCond,x,t)
    exactSolution1b = lambda x,t: exactSolution1a(x,t)*exp(-mu*4.*t*pi**2)

    sol = None
    if type == 'a':
        sol = exactSolution1a
    elif type == 'b':
        sol = exactSolution1b

    ssol = lambda x: sol(x,T)
    u0 = my_array_copy(my_vector_vectorize(initCond,X))
    exact = my_array_copy(my_vector_vectorize(ssol,X))
    # u0 = np.array([initCond(x) for x in X],dtype=np.float64)
    # exact = np.array([sol(x,T) for x in X],dtype=np.float64)

    u = solver(u0, M)
    return u, exact, X

# You can insert a different letters (a or b) to get different plots.
# TEST 1
for scheme, name in zip([upwind, laxWendroff],['upwind', 'LaxWendroff']):
    grids = [10**i for i in range(1, 4)]
    errors = []
    for n in grids:
        t = problemSetup(scheme,initialCondition1,'b',n,nu,a)
        errors.append(maxNorm(t[0],t[1]))
    plt.loglog(grids, errors, label=name)
plt.legend(loc='best')
plt.show()

# Comment the previous section and uncomment this. Run again and see three different plots for different grids
# TEST 2
# grids = [10**i for i in range(1, 4)]
# for n in grids:
#     plt.figure()
#     plt.title(str(n))
#     for scheme, name in zip([upwind, laxWendroff],['upwind', 'LaxWendroff']):
#         t = problemSetup(scheme,initialCondition2,'a',n,nu,a)
#         plt.plot(t[2],t[0],label=name)
#     plt.plot(t[2],t[1],label='exact')
#     plt.margins(x=0)
#     plt.legend(loc='best')
# plt.show()
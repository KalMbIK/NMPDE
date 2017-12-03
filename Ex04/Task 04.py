import numpy as np
import sys
import math
from numpy import sin, pi, fabs, exp
import matplotlib.pyplot as plt

def init_cond1(_X):
    return sin(2. * pi * _X)

def exactA(_X, _t):
    return sin(2.*pi*(_X-_t))

def exactB(_X, _t, _h, nu):
    mu=_h/2.*(1.-fabs(nu))
    return sin(2.*pi*(_X-_t))*exp(-mu*4.*pi*pi*_t)

def error(numeric, exact):
    if len(numeric)!=len(exact):
        return 'Incompatible sizes!'
    max=np.fabs(numeric[0]-exact[0])
    for i in range(1, len(exact)):
        temp=np.fabs(numeric[i]-exact[i])
        if temp>max:
            max=temp
    return max

def scheme(init, _X, _T, nu, flag):
    if flag=='uw':
        s=fabs(nu)
    elif flag == 'lw':
        s=nu*nu
    else:
        sys.exit('Wrong flag!')
    u=init
    dim=len(u)
    for n in range(len(_T)-1):
        u1=np.zeros(dim)
        u1[0]=u[0]-nu/2.*(u[1]-u[dim-1])+s/2.*(u[1]-2.*u[0]+u[dim-1])
        for i in range(1, dim-1):
            u1[i]=u[i]-nu/2.*(u[i+1]-u[i-1])+s/2.*(u[i+1]-2.*u[i]+u[i-1])
        u1[dim-1]=u[dim-1]-nu/2.*(u[0]-u[dim-2])+s/2.*(u[0]-2.*u[dim-1]+u[dim-2])
        u=u1
    return u

N=100
nu=0.9
h=1./N
tau=nu*h
M=int(N/nu+1)
Tc=M*tau
X = np.arange(0., 1., h)
T = np.arange(0., Tc + tau, tau)

# num=scheme(init_cond1(X), X, T, nu, 'lw')
# exa=exactB(X, Tc, h, nu)

# plt.plot(X, num)
# plt.plot(X, exa)
# plt.show()
num=scheme(init_cond1(X), X, T, nu, 'uw')
print 'Error in 1A for N='+str(N)+': '+str(error(num, exactA(X, Tc)))
print 'Error in 1B for N='+str(N)+': '+str(error(num, exactB(X, Tc, h, nu)))

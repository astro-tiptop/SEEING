import numpy as np
import cupy as cp

import scipy
import scipy.special

import cupyx
import cupyx.scipy
import cupyx.scipy.special
import cupyx.scipy.linalg


cpulibExt = {
    'factorial': scipy.special.factorial,
    'binomial': scipy.special.binom,
    'besselj': scipy.special.jv,
    'besselk': scipy.special.kv,
    'besseli': scipy.special.iv,
    'bessely': scipy.special.yv,
    'erf': scipy.special.erf,
    'gamma': scipy.special.gamma
}

cpulib = ["scipy", cpulibExt] 

#
# bessel j per ordine non intero, da fare
#

def besselj__n(n, z):
    if n<0:
        return -1**(-n) * besselj__n(-n, z)
    if n==0:
        return cupyx.scipy.special.j0(z)
    elif n==1:
        return cupyx.scipy.special.j1(z)
    elif n>=2:
        return 2*besselj__n(int(n)-1, z)/z - besselj__n(int(n)-2, z)


def besselk__n(a, z):
    return (besselj__n(a,z)*cp.cos(a*cp.pi) - besselj__n(-a,z) ) / cp.sin(a*cp.pi)
    
def gamma__1(n):
    return cupyx.scipy.special.gamma(n+1)


gpulib = {
    'pow': cp.power,
    'atan': cp.arctan,
    'atan2': cp.arctan2,
    'I': 1j,
    'pi': cp.pi,
    'besselj': besselj__n,
    'besselk': besselk__n,
    'factorial': gamma__1,
    'gamma': cupyx.scipy.special.gamma,
    'erf': cupyx.scipy.special.erf
}


for key, value in cp.__dict__.items():
    gpulib[key] = value

backendLibs = {
    cp: gpulib,
    np: cpulib
}

'''
Not necessary: 

    'sqrt': cp.sqrt,
    'cos': cp.cos,
    'sin': cp.sin,
    'exp': cp.exp,
'''
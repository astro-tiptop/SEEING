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


def gamma__1(n):
    return cupyx.scipy.special.gamma(n+1)


factorial__ = gamma__1

# direct and inverse trigonometric functions


def cot__(z):
    return 1.0/cp.tan(z)


def sec__(z):
    return 1.0/cp.cos(z)


def csc__(z):
    return 1.0/cp.sin(z)


def coth__(z):
    return 1.0/cp.tanh(z)


def arccot__(z):
    return cp.arctan(1.0/z)


def arcsec__(z):
    return cp.arccos(1.0/z)


def arccsc__(z):
    return cp.arcsin(1.0/z)


# direct and inverse hyperbolic functions


def sech__(z):
    return 1.0/cp.cosh(z)


def csch__(z):
    return 1.0/cp.sinh(z)


def arccoth__(z):
    return cp.arctanh(1.0/z)


def arcsech__(z):
    return cp.arccosh(1.0/z)


def arccsch__(z):
    return cp.arcsinh(1.0/z)


def frac__(z):
    return z+1-cp.ceil(z)


def identity__():
    return z


def root__(z, n):
    return z**(1.0/n)


def cbrt__(z):
    return z**(1.0/3.0)


def binomial_coeff__(n, k):
    return factorial__(n)/(factorial__(n-k)*factorial__(k))


def bell__(n):
    if n==0:
        return 1
    else:
        s=0
        for k in range(n):
            s+= bell__(k) * binomial_coeff__(n-1, k)
        return s

    
def catalan__(n):
    return binomial_coeff__(2*n, n)/(n+1)


def subfactorial__(n):
    return cp.ceil(cp.round_(factorial__(n)/cp.exp(1))) - 1.0


def falling_factorial__(x, k):
    return ((-1)**k)*cupyx.scipy.special.gamma(k - x)/cupyx.scipy.special.gamma(-x)


def fibonacci__(n):
    if n==0:
        return 0
    elif n==1:
        return 1
    else:
        return fibonacci__(n-1) + fibonacci__(n-2)

    
def tribonacci__(n):
    if n==0:
        return 0
    elif n==1:
        return 1
    elif n==2:
        return 1
    else:
        return tribonacci__(n-1) + tribonacci__(n-2) + tribonacci__(n-3)

    
def lucas__(n):
    if n==0:
        return 2
    elif n==1:
        return 1
    else:
        return lucas__(n-1) + lucas__(n-2)

    
def rising_factorial__(n,k):
    return factorial__(k + n - 1)/factorial__(n - 1)


def digamma__(z):
    return cupyx.scipy.special.polygamma(0, z)


def trigamma__(z):
    return cupyx.scipy.special.polygamma(1, z)


def besselj__n(n, z):
    if n<0:
        return -1**(-n) * besselj__n(-n, z)
    if n==0:
        return cupyx.scipy.special.j0(cp.real(z))
    elif n==1:
        return cupyx.scipy.special.j1(cp.real(z))
    elif n>=2:
        return 2*(n-1)*besselj__n(int(n)-1, z)/cp.real(z) - besselj__n(int(n)-2, z)

    
def bessely__n(n, z):
    return (besselj__n(n,z)*cp.cos(n*cp.pi) - besselj__n(-n,z) ) / cp.sin(n*cp.pi)

    
def besselk__n(n, z):
    return cp.pi/2 * (besseli__n(n, z)-besseli__n(-n, z))/cp.sin(n*cp.pi)
    
    
def besseli__n(n, z):
    return (1j**-n) * besselj__n(n, 1j*z)


def hankel1__n(n, z):
    return besselj__n(n, z) + 1j * bessely__n(n, z)


def hankel2__n(n, z):
    return besselj__n(n, z) - 1j * bessely__n(n, z)

def floatArray(z):
    return cp.array(z, cp.float64)

gpulib = {
    'ImmutableDenseMatrix': floatArray,
    'E': np.exp(1.0),
    'I': 1j,
    'pi': cp.pi,
    're': cp.real,
    'im': cp.imag,
    'arg': cp.angle,
    'conjugate': cp.conj,
#    'polar_lift': cp. ,
#    'periodic_argument': cp.,
#    'principal_branch': cp.,
    'asin': cp.arcsin,
    'acos': cp.arccos,
    'atan': cp.arctan,
    'asinh': cp.arcsinh,
    'acosh': cp.arccosh,
    'atanh': cp.arctanh,
    'cot': cot__,
    'sec': sec__,
    'csc': csc__,
    'acot': arccot__,
    'asec': arcsec__,
    'acsc': arccsc__,
    'atan2': cp.arctan2,
    'coth': coth__,
    'sech': sech__,
    'csch': csch__,
    'acoth': arccoth__,
    'asech': arcsech__,
    'acsch': arccsch__,
    'ceiling': cp.ceil,
    'RoundFunction': cp.round_,
    'frac': frac__,
#    'LambertW': ,
#    'exp_polar': ,
#    'ExprCondPair': ,
#    'Piecewise': cp.piecewise, # ? unlikely working
    'Identity': identity__,
    'Min': cp.fmin,
    'Max': cp.fmax,
    'root': root__, 
    'cbrt': cbrt__,
#    'real_root': ,
    'bell': bell__,
#    'bernoulli': bernoulli__,
    'binomial': binomial_coeff__,
    'catalan': catalan__,
 #   'euler': euler__,
    'factorial': factorial__,
    'subfactorial': subfactorial__,
#    'factorial2': ,
    'FallingFactorial': falling_factorial__,
    'fibonacci': fibonacci__,
    'tribonacci': tribonacci__,
#    'harmonic': ,
    'lucas': lucas__,
#    'genocchi': ,
#    'partition': partition__,
#    'MultiFactorial': ,
    'RisingFactorial': rising_factorial__,
#    'stirling': ,
#    'nC': ,
#    'nP': ,
#    'nT': ,    
#    'Heaviside': heaviside__,
    'gamma': cupyx.scipy.special.gamma,
    'loggamma': cupyx.scipy.special.gammaln,
    'polygamma': cupyx.scipy.special.polygamma,
    'digamma': digamma__,
    'trigamma': trigamma__,
#    'uppergamma': uppergamma__,
#    'lowergamma': lowergamma__,
#    'multigamma': lowergamma__,
    'erf': cupyx.scipy.special.erf,
    'erfc': cupyx.scipy.special.erfc,
#    'erfi': erfi__,
#    'erf2': erf2__,
    'erfinv': cupyx.scipy.special.erfinv,
    'erfcinv': cupyx.scipy.special.erfcinv,
#    'erf2inv': erf2inv__,
    'pow': cp.power,
    'besselj': besselj__n,
    'bessely': bessely__n,
    'besseli': besseli__n,
    'besselk': besselk__n,
    'hankel1': hankel1__n,
    'hankel2': hankel2__n,    
    'zeta': cupyx.scipy.special.zeta
}


# this will add a lot of junk to the dictionary
# but also all the obvios translations, like 'sin':(cp.)sin
# the junk wil not be a problem, since the relative sympy functions do not exists
# items like: 'atleast_1d':atleast_1d 
for key, value in cp.__dict__.items():
    gpulib[key] = value


backendLibs = {
    cp: gpulib,
    np: cpulib
}

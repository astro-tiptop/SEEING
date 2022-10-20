import scipy
import scipy.special
import numpy as np

from . import gpuEnabled

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

cp = np
gpulib = cpulib

if gpuEnabled:
    from .lambdifyDictionariesGPU import *

backendLibs = {
    cp: gpulib,
    np: cpulib
}

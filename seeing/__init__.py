import numpy as np

gpuEnabled = False
cp = None

try:
    import cupy as cp
    print("Cupy import successfull. Installed version is:", cp.__version__)
    gpuEnabled = True
except:
    print("Cupy import failed. SEEING will fall back to CPU use.")
    cp = np

from seeing.sympyHelpers import *
from seeing.formulary import *
from seeing.integrator import *
from seeing._version import __version__

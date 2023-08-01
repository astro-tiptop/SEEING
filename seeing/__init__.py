import numpy as np
import os

gpuEnabled = False
cp = None

systemDisable = os.environ.get('SEEING_DISABLE_GPU', 'FALSE')
if systemDisable=='FALSE':
    try:
        import cupy as cp
        print("Cupy import successfull. Installed version is:", cp.__version__)
        gpuEnabled = True
    except:
        print("Cupy import failed. SEEING will fall back to CPU use.")
        cp = np
else:
    print("env variable SEEING_DISABLE_GPU prevents using the GPU.")
    cp = np

from seeing.sympyHelpers import *
from seeing.formulary import *
from seeing.integrator import *
from seeing._version import __version__

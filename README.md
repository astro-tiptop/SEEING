# SEEING

Sympy Expressions Evaluation Implemented oN the GPU

The idea is to provide the following to Sympy users:

- A simple way to evaluate an expressions on the GPU
- Backend agnostic evaluation (currently numpy and cupy are supported, in the future clpy or other numerical backends might be added)
- Support tools to handle groups of expressions (evaluate, plot, substitute parameters, share variables)
- Numerical methods (for now integration over n-dimensional domains)

Numpy/Scipy evaluation comes out of the box with Sympy, using the lambdify
feature.  In order to do the same on the GPU, first development effort of
SEEING relies Cupy only, thus you need to have a CUDA enabled GPU to use it.

## Installation

SEEING can be installed with pip:

    pip install astro-seeing

For GPU usage SEEING relies on Cupy, which can be installed as a dependency
with:

    pip install "astro-seeing[gpu]"

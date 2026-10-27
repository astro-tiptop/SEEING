import numpy as np
import sympy as sp
import pytest

from seeing import cp, gpuEnabled
from seeing.integrator import Integrator, _cached_lambdify
from seeing.lambdifyDictionaries import backendLibs

A_VALUES = np.array([0.5, 1.0, 2.0, 4.0])
L = 10.0
N_POINTS = 2001


def _gaussian_problem():
    """Build F(a) = Integral(exp(-a x^2), (x, -L, L)) and its sampling spec."""
    x, a = sp.symbols('x a', real=True)
    lh = sp.Function('F')(a)
    integral = sp.Integral(sp.exp(-a * x**2), (x, -L, L))
    params = [('a', A_VALUES, A_VALUES, N_POINTS, 'provided')]
    return lh, integral, params


def _new_integrator(xp=np):
    # 'plain' returns the raw integral ('intensity' would square it)
    return Integrator(xp, float, integrationMode='plain')


requires_gpu = pytest.mark.skipif(not gpuEnabled, reason='cupy/GPU not available')


@pytest.mark.parametrize('method', ['trap', 'rect'])
def test_gaussian_integral(method):
    lh, integral, params = _gaussian_problem()
    _, res = _new_integrator().IntegralEval(lh, integral, params, method=method)
    np.testing.assert_allclose(np.ravel(res), np.sqrt(np.pi / A_VALUES), rtol=1e-4)


def test_gaussian_integral_equation_form():
    lh, integral, params = _gaussian_problem()
    _, res = _new_integrator().IntegralEvalE(sp.Eq(lh, integral), params, method='trap')
    np.testing.assert_allclose(np.ravel(res), np.sqrt(np.pi / A_VALUES), rtol=1e-4)


def test_cache_hits_across_calls_and_instances():
    _cached_lambdify.cache_clear()
    lh, integral, params = _gaussian_problem()
    _, r1 = _new_integrator().IntegralEval(lh, integral, params, method='trap')
    info1 = _cached_lambdify.cache_info()
    assert info1.misses == 1
    # same instance type, new instance, same backend
    _, r2 = _new_integrator().IntegralEval(lh, integral, params, method='trap')
    info2 = _cached_lambdify.cache_info()
    assert info2.hits > info1.hits
    assert info2.misses == info1.misses
    assert np.array_equal(r1, r2)


def test_function_eval_cache():
    _cached_lambdify.cache_clear()
    x = sp.symbols('x')
    f = sp.sin(x) ** 2
    params = [('x', 0.0, 3.0, 50, 'linear')]
    _, y1 = _new_integrator().functionEval(f, params)
    info1 = _cached_lambdify.cache_info()
    _, y2 = _new_integrator().functionEval(f, params)
    info2 = _cached_lambdify.cache_info()
    assert info2.hits == info1.hits + 1
    assert info2.misses == info1.misses
    assert np.array_equal(y1, y2)
    xs = np.linspace(0.0, 3.0, 50)
    np.testing.assert_allclose(np.ravel(y1), np.sin(xs) ** 2, rtol=1e-12)


def test_cached_lambdify_matches_sympy():
    x, y = sp.symbols('x y')
    expr = sp.besselj(1, x) * y + sp.exp(-y)
    rng = np.random.default_rng(0)
    xv = rng.uniform(0.1, 5.0, 100)
    yv = rng.uniform(0.1, 5.0, 100)
    ref = sp.lambdify((x, y), expr, backendLibs[np])
    cached = _cached_lambdify((x, y), expr, np)
    assert np.array_equal(cached(xv, yv), ref(xv, yv))


def test_list_args_rejected():
    # args are part of the lru_cache key: callers must pass a tuple
    x, y = sp.symbols('x y')
    with pytest.raises(TypeError):
        _cached_lambdify([x, y], x + y, np)


def test_no_cross_contamination():
    x, y = sp.symbols('x y')
    xv, yv = np.array([3.0, 5.0]), np.array([1.0, 7.0])
    f_xy = _cached_lambdify((x, y), x - y, np)
    f_yx = _cached_lambdify((y, x), x - y, np)
    assert f_xy is not f_yx
    assert np.array_equal(f_xy(xv, yv), xv - yv)
    assert np.array_equal(f_yx(xv, yv), -(xv - yv))
    # different expression, same args
    f_sum = _cached_lambdify((x, y), x + y, np)
    assert np.array_equal(f_sum(xv, yv), xv + yv)
    # earlier entries are unchanged
    assert np.array_equal(_cached_lambdify((x, y), x - y, np)(xv, yv), xv - yv)


@requires_gpu
def test_backend_isolation_in_cache():
    _cached_lambdify.cache_clear()
    x = sp.symbols('x')
    expr = sp.besselj(0, x) * sp.exp(-x)
    xs = np.linspace(0.1, 3, 50)
    f_cp = _cached_lambdify((x,), expr, cp)
    f_np = _cached_lambdify((x,), expr, np)
    info = _cached_lambdify.cache_info()
    assert info.misses == 2 and info.hits == 0
    assert f_cp is not f_np
    # repeated lookups hit their own entry
    assert _cached_lambdify((x,), expr, cp) is f_cp
    assert _cached_lambdify((x,), expr, np) is f_np
    y_cp = f_cp(cp.asarray(xs))
    y_np = f_np(xs)
    assert isinstance(y_cp, cp.ndarray)
    assert isinstance(y_np, np.ndarray)
    np.testing.assert_allclose(cp.asnumpy(y_cp), y_np, rtol=1e-12, atol=1e-14)


@requires_gpu
@pytest.mark.parametrize('order', ['numpy_first', 'cupy_first'])
def test_numpy_and_cupy_integrators_coexist(order):
    _cached_lambdify.cache_clear()
    lh, integral, params = _gaussian_problem()
    backends = [np, cp] if order == 'numpy_first' else [cp, np]
    integrators = [_new_integrator(xp) for xp in backends]
    # two rounds so that the second one is served from the cache
    for _ in range(2):
        for xp, integrator in zip(backends, integrators):
            _, res = integrator.IntegralEval(lh, integral, params, method='trap')
            # IntegralEval returns host arrays for both backends
            np.testing.assert_allclose(np.ravel(cp.asnumpy(res)),
                                       np.sqrt(np.pi / A_VALUES), rtol=1e-4)
    assert _cached_lambdify.cache_info().misses == 2

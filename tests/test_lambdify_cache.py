import numpy as np
import sympy as sp
import pytest

from seeing.integrator import Integrator, _cached_lambdify, _lambdify_by_id
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


def _new_integrator():
    # 'plain' returns the raw integral ('intensity' would square it)
    return Integrator(np, float, integrationMode='plain')


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
    _lambdify_by_id.cache_clear()
    lh, integral, params = _gaussian_problem()
    _, r1 = _new_integrator().IntegralEval(lh, integral, params, method='trap')
    info1 = _lambdify_by_id.cache_info()
    assert info1.misses == 1
    # same instance type, new instance, same backend
    _, r2 = _new_integrator().IntegralEval(lh, integral, params, method='trap')
    info2 = _lambdify_by_id.cache_info()
    assert info2.hits > info1.hits
    assert info2.misses == info1.misses
    assert np.array_equal(r1, r2)


def test_function_eval_cache():
    _lambdify_by_id.cache_clear()
    x = sp.symbols('x')
    f = sp.sin(x) ** 2
    params = [('x', 0.0, 3.0, 50, 'linear')]
    _, y1 = _new_integrator().functionEval(f, params)
    info1 = _lambdify_by_id.cache_info()
    _, y2 = _new_integrator().functionEval(f, params)
    info2 = _lambdify_by_id.cache_info()
    assert info2.hits == info1.hits + 1
    assert info2.misses == info1.misses
    assert np.array_equal(y1, y2)
    xs = np.linspace(0.0, 3.0, 50)
    np.testing.assert_allclose(np.ravel(y1), np.sin(xs) ** 2, rtol=1e-12)


@pytest.mark.parametrize('as_type', [list, tuple])
def test_cached_lambdify_matches_sympy(as_type):
    x, y = sp.symbols('x y')
    expr = sp.besselj(1, x) * y + sp.exp(-y)
    modules = backendLibs[np]
    rng = np.random.default_rng(0)
    xv = rng.uniform(0.1, 5.0, 100)
    yv = rng.uniform(0.1, 5.0, 100)
    ref = sp.lambdify((x, y), expr, modules)
    cached = _cached_lambdify(as_type([x, y]), expr, modules)
    assert np.array_equal(cached(xv, yv), ref(xv, yv))


def test_list_and_tuple_args_share_entry():
    _lambdify_by_id.cache_clear()
    x, y = sp.symbols('x y')
    expr = x + y
    modules = backendLibs[np]
    f1 = _cached_lambdify([x, y], expr, modules)
    f2 = _cached_lambdify((x, y), expr, modules)
    assert f1 is f2
    assert _lambdify_by_id.cache_info().misses == 1


def test_no_cross_contamination():
    x, y = sp.symbols('x y')
    modules = backendLibs[np]
    xv, yv = np.array([3.0, 5.0]), np.array([1.0, 7.0])
    f_xy = _cached_lambdify((x, y), x - y, modules)
    f_yx = _cached_lambdify((y, x), x - y, modules)
    assert f_xy is not f_yx
    assert np.array_equal(f_xy(xv, yv), xv - yv)
    assert np.array_equal(f_yx(xv, yv), -(xv - yv))
    # different expression, same args
    f_sum = _cached_lambdify((x, y), x + y, modules)
    assert np.array_equal(f_sum(xv, yv), xv + yv)
    # earlier entries are unchanged
    assert np.array_equal(_cached_lambdify((x, y), x - y, modules)(xv, yv), xv - yv)

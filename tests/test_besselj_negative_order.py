import numpy as np
import pytest
from scipy.special import jv

cp = pytest.importorskip('cupy')

try:
    cp.cuda.runtime.getDeviceCount()
except Exception:
    pytest.skip('no CUDA device', allow_module_level=True)

from seeing.lambdifyDictionariesGPU import besselj__n


@pytest.mark.parametrize('n', range(-3, 4))
def test_besselj_integer_order_matches_scipy(n):
    z = np.linspace(0.1, 20.0, 400)
    out = cp.asnumpy(besselj__n(n, cp.asarray(z)))
    np.testing.assert_allclose(out, jv(n, z), rtol=0, atol=1e-12)


def test_besselj_negative_order_parity():
    z = cp.linspace(0.5, 20.0, 100)
    for n in (1, 2, 3):
        np.testing.assert_allclose(cp.asnumpy(besselj__n(-n, z)),
                                   (-1)**n * cp.asnumpy(besselj__n(n, z)), rtol=0, atol=0)

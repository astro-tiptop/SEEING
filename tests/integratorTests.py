DEV_PATH = '/home/frossi/dev/'
import sys
if not DEV_PATH in sys.path:
    sys.path.append(DEV_PATH)

    
from SEEING.seeing import *
import unittest


class TestEvaluateFormula(unittest.TestCase):
    def test_simple_expression(self):
        """
        Test 
        """
        aname1 = 'x'
        aname2 = 'y'
        x = sp.symbols(aname1)
        y = sp.symbols(aname2)
        aexpr = x**2 + y**2
        result1 = evaluateFormula( aexpr, {}, ['x', 'y'], [2.0*np.ones(10), 3.0*np.ones(20)] )
        self.assertTrue(np.allclose(result1[2], 13.0*np.ones((20,10))))
        result2 = evaluateFormula( aexpr, {}, ['x', 'y'], [2.0*np.ones((20,10)), 3.0*np.ones((20,10))] )
        self.assertTrue(np.allclose(result2[2], 13.0*np.ones((20,10))))


'''
x, y = sp.symbols('x y')
afunc = 1 - x**2 - y**2
paramsAndRange = [( 'x', -1, 1, 1000, 'linear' ), ( 'y', -1, 1, 1000, 'linear' )]
vplot, zplot = mIt.functionEval(afunc, paramsAndRange)
plt.imshow(zplot)
'''


        
if __name__ == '__main__':
    unittest.main()
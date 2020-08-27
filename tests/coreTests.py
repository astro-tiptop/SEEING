DEV_PATH = '/home/frossi/dev/'
import sys
if not DEV_PATH in sys.path:
    sys.path.append(DEV_PATH)

from SEEING.seeing import *
import unittest


class TestGetSymbolByName(unittest.TestCase):
    def test_simple_expression(self):
        """
        Test 
        """
        aname = 'x'
        x = sp.symbols(aname)
        aexpr = x ** 2
        result = getSymbolByName(aexpr, aname)
        self.assertEqual(result, x)

        
class TestSubsParamsByName(unittest.TestCase):
    def test_simple_expression(self):
        """
        Test 
        """
        aname1 = 'x'
        aname2 = 'y'
        x = sp.symbols(aname1)
        y = sp.symbols(aname2)
        aexpr = x ** 2
        result = subsParamsByName(aexpr, {aname1:y})
        self.assertEqual(getSymbolByName(result, aname2), y)

        
class TestGetRestrictedLambdaBasic(unittest.TestCase):
    def test_simple_expression(self):
        """
        Test 
        """
        aname1 = 'x'
        aname2 = 'y'
        x = sp.symbols(aname1)
        y = sp.symbols(aname2)
        aexpr = x**2 + y**2
        alambda = getRestrictedLambdaBasic(aexpr, {aname1:1}, ['y'])
        result = alambda(2*np.ones(10))
        self.assertTrue(np.allclose(result, 5*np.ones(10)))

        
class TestEvaluateLambda(unittest.TestCase):
    def test_simple_expression(self):
        """
        Test 
        """
        aname1 = 'x'
        aname2 = 'y'
        x = sp.symbols(aname1)
        y = sp.symbols(aname2)
        aexpr = x**2 + y**2
        alambda = getRestrictedLambdaBasic(aexpr, {}, ['x', 'y'])
        result1 = evaluateLambda( alambda, ['x', 'y'], [2.0*np.ones(10), 3.0*np.ones(20)] )
        self.assertTrue(np.allclose(result1[2], 13.0*np.ones((20,10))))
        result2 = evaluateLambda( alambda, ['x', 'y'], [2.0*np.ones((20,10)), 3.0*np.ones((20,10))] )
        self.assertTrue(np.allclose(result2[2], 13.0*np.ones((20,10))))

        
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

        
        
        
        
if __name__ == '__main__':
    unittest.main()
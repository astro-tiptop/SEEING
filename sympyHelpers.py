import scipy.ndimage
import scipy.interpolate
import sympy as sp
import inspect

import matplotlib.pyplot as plt

import numpy as np
import scipy
import scipy.signal

import cupy as cp
import cupyx
import cupyx.scipy
import cupyx.scipy.special

from functools import reduce

#mempool = cp.get_default_memory_pool()
# print(mempool.used_bytes())              # 0
# print(mempool.total_bytes())             # 0

pplib = {
    'factorial': scipy.special.factorial,
    'binomial': scipy.special.binom,
    'besselj': scipy.special.jv,
    'besselk': scipy.special.kv,
    'besseli': scipy.special.iv,
    'bessely': scipy.special.yv,
    'erf': scipy.special.erf,
    'gamma': scipy.special.gamma
}

nplib = ["scipy", pplib] 

def bessel__n(n, z):
    if n==0:
        return cupyx.scipy.special.j0(z)
    elif n==1:
        return cupyx.scipy.special.j1(z)
    elif n>=2:
        return 2*bessel__n(n-1, z)/z - bessel__n(n-2, z)

def gamma__1(n):
    return cupyx.scipy.special.gamma(n+1)
    
cplib = {
    'sqrt': cp.sqrt,
    'cos': cp.cos,
    'sin': cp.sin,
    'exp': cp.exp,
    'pow': cp.power,
    'atan': cp.arctan,
    'atan2': cp.arctan2,
    'I': 1j,
    'pi': cp.pi,
    'besselj': bessel__n,
    'factorial': gamma__1,
    'gamma': cupyx.scipy.special.gamma,
    'erf': cupyx.scipy.special.erf
}

for key, value in cp.__dict__.items():
    cplib[key] = value


def getSymbolByName(expr, sname):
    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
    return free_syms[sname]


def subsParamsByName(expr1, subsDict):
    expr = expr1
#    syms_and_funs = set(expr.free_symbols) | set([i.func for i in expr.atoms(sp.Function) if isinstance(i, sp.function.AppliedUndef)])
#    free_syms = {symbol.name: symbol for symbol in syms_and_funs}
    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
#    print(free_syms)
    for sname, sval in subsDict.items():
        expr = expr.subs(free_syms[sname], sval)
    return expr


def replaceParamsByName(expr, subsDict):    
    syms_and_funs = set(expr.free_symbols) | set([i.func for i in expr.atoms(sp.Function) if isinstance(i, sp.function.AppliedUndef)])
    free_syms = {symbol.name: symbol for symbol in syms_and_funs}
#    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
    for sname, sval in subsDict.items():
#        print(sname, sval)
#        print(type(free_syms[sname]))
        expr = expr.replace({free_syms[sname], sval})
    return expr


def _getRestrictedLambda(
        expr,
        subsDict,
        independentVarNames,
        lmodules='scipy'):
    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
    for sname, sval in subsDict.items():
        expr = expr.subs(free_syms[sname], sval)
    ll = []
    ii = 0
    for vname in independentVarNames:
        if vname in free_syms.keys():
            ll.append(free_syms[vname])
        else:
            ll.append(sp.Symbol('dummy' + str(ii)))
            ii += 1
    tt = tuple(ll)
    return sp.lambdify(tt, expr, modules=lmodules)


def evaluateLambda(p_f, independentVarNames, samples):

    #print(inspect.getsource(p_f))

    nvars = len(independentVarNames)
    if nvars == 1:
        return (samples[0], p_f(samples[0]))
    elif nvars == 2 and samples[0].ndim == 1 and samples[1].ndim == 1:
        (v1, v2) = np.meshgrid(samples[0], samples[1])
        return (v1, v2, p_f(v1, v2))
    elif nvars == 2 and samples[0].ndim == 2 and samples[1].ndim == 2:
        return (samples[0], samples[1], p_f(samples[0], samples[1]))
    else:
        return None


def evaluateFormula(
        expr,
        subsDict,
        independentVarNames,
        samples,
        modules='scipy'):
    p_f = _getRestrictedLambda(expr, subsDict, independentVarNames, modules)
#    print(inspect.getsource(p_f))
    return evaluateLambda(p_f, independentVarNames, samples)


class Formulary(object):
    # each element of the formulas list is a tuple (_lhs, _rhs, _eq)
    def __init__(self, name='', names=[], formulas=[]):
        self.name = name
#        self.variables = variables
        self.formulas = dict(zip(names, formulas))

    def addFormula(self, name, formula):
        self.formulas[name] = formula

    def displayAll(self):
        for name, f in self.formulas.items():
            #            print(name)
            display(f[2])

    def getFormula(self, name):
        return self.formulas[name]

    def getFormulaRhs(self, name):
        return self.formulas[name][1]

    def getFormulaLhs(self, name):
        return self.formulas[name][0]

    def getRestrictedLambda(
            self,
            name,
            subsDict,
            independentVarNames,
            modules='scipy'):
        _, _expr, _ = self.getFormula(name)
        return _getRestrictedLambda(
            _expr, subsDict, independentVarNames, modules)

    def evaluateFormula(
            self,
            name,
            subsDict,
            independentVarNames,
            samples,
            modules='scipy'):
        p_f = self.getRestrictedLambda(
            name, subsDict, independentVarNames, modules)
        return evaluateLambda(p_f, independentVarNames, samples)

    def plotFormula(
            self,
            name,
            subsDict,
            independentVarNames,
            samples,
            modules='scipy',
            log=False):
        _data = self.evaluateFormula(
            name, subsDict, independentVarNames, samples, modules)
        nvars = len(independentVarNames)
        if nvars == 1 and samples[0].ndim == 1:
            fig = plt.figure(figsize=(7, 7))
            axs = fig.add_subplot(111)
            axs.axis('auto')
            if log:
                plt.yscale('log')
                plt.xscale('log')
            else:
                plt.yscale('linear')
                plt.xscale('linear')
            axs.plot(*_data)
            #print(_data)
        elif nvars == 1 and samples[0].ndim == 2:
            plt.imshow(_data[1], cmap='jet')
        elif nvars == 2:
            plt.imshow(_data[2], cmap='jet')
        else:
            pass
        plt.show()
        

class Integrator(object):
    # mode can be:
    # plain: simply evaluate the integral
    # absolute: perfeorms also absolute value ( modulus for complex integrals) at each evaluation point
    # intensity: perfeorms also square of absolute value ( modulus for complex
    # integrals) at each evaluation point
    def __init__(
            self,
            mode='intensity',
            xp=np,
            evalType=float):
        self.xp = xp
        if self.xp==cp:
            self.modules = cplib
        else:
            self.modules = nplib
        
        self.evalType = evalType
        if mode == 'intensity':
            self.intensityMap = self.xp.square
            self.postMap = self.xp.absolute
        elif mode == 'absolute':
            self.intensityMap = lambda x: x
            self.postMap = self.xp.absolute
        else:
            self.intensityMap = lambda x: x
            self.postMap = lambda x: x

    def getSampling(self, l, h, npoints, spacing='linear'):
        if spacing == 'geometric':
            return self.xp.asarray(np.geomspace(l, h, npoints), dtype=self.evalType)
        elif spacing == 'sqrt':
            return self.xp.asarray(h*np.sqrt(np.linspace(l/h, h/h, npoints), dtype=self.evalType))
        elif spacing == 'random':
            return l + self.xp.random.random(npoints, dtype=self.evalType)*(h-l)
        else:
            return self.xp.asarray(np.linspace(l, h, npoints), dtype=self.evalType)

    def outputData(self, _data):
        if self.xp == cp:
            return cp.asnumpy(_data)
        else:
            return _data
    def getWeightsArray(self, n):
        p0 = np.ones(n, dtype=np.float64)
        p0[0] = 0.5
        p0[-1] = 0.5
        return np.asarray(p0, dtype=np.float64)        
        
    def parametricIntegralEvaluation(self, integrationVarsSamplings, paramsSamplings, integrandFunction, method='trap', smallSamplings=[]):
        integrationAndParamsVarSamplingGrids = self.xp.meshgrid( *[*integrationVarsSamplings, *paramsSamplings], sparse=True, copy=False)
        weightsArrays = [self.getWeightsArray(aa.shape[0]) for aa in integrationVarsSamplings]
        weightsGrid = reduce(np.multiply.outer, weightsArrays[::-1])
        weightsGrid = self.xp.asarray( np.dstack([weightsGrid] * paramsSamplings[0].shape[0]) )
        
        nVars = len(integrationVarsSamplings)
        nParams = len(paramsSamplings)
        # this is strange, but for some reason tt=(1)  
        # instead of tt=(0) when we have one integration variable
        # while when we have more than one tt=(0,1,..,n)
        tt = (1)
        if nVars>1:
            tt=tuple(range(nVars))

        def genericSum(xx):
            return self.xp.sum( xx, axis=tt )
            
        def genericReduction(xx):
            return self.postMap( self.xp.sum( xx * weightsGrid, axis=tt ) )
        
        if self.xp == cp:
            if method=='rect':
                @cp.fuse(kernel_name='integratedFunctionRect')
                def integratedFunction(*integrationAndParamsVarSamplingGrids, reduce=genericSum, post_map=self.postMap):                               
                    return post_map(reduce( integrandFunction(*integrationAndParamsVarSamplingGrids))) 
            else:
                @cp.fuse(kernel_name='integratedFunctionTrap')
                def integratedFunction(*integrationAndParamsVarSamplingGrids):                               
                    return integrandFunction(*integrationAndParamsVarSamplingGrids) 
        else:
            integrandFunctionV = np.vectorize(integrandFunction)
            if method=='rect':
                def integratedFunction(*integrationAndParamsVarSamplingGrids, post_map=self.postMap):
    #                return post_map(genericSum(np.nan_to_num(integrandFunction(*integrationAndParamsVarSamplingGrids))))
                    return post_map(genericSum(integrandFunctionV(*integrationAndParamsVarSamplingGrids)))
            else:
                def integratedFunction(*integrationAndParamsVarSamplingGrids):
                    return integrandFunctionV(*integrationAndParamsVarSamplingGrids)
            
        scaleFactor = 1.0
        for ii in range(nVars):
            if (len(integrationVarsSamplings[ii])>1):
                scaleFactor *= (integrationVarsSamplings[ii][1] - integrationVarsSamplings[ii][0])
        if method=='trap':
            return scaleFactor * genericReduction(integratedFunction(*integrationAndParamsVarSamplingGrids))
        elif method=='rect':
            return scaleFactor * integratedFunction(*integrationAndParamsVarSamplingGrids)
#        elif method=='mc':
#            integrationAndParamsVarSamplingGridsSmall = self.xp.meshgrid( *[*smallSamplings, *paramsSamplings], sparse=True, copy=False)
#            evalf0 = integratedFunction(*integrationAndParamsVarSamplingGridsSmall)
#            f_max = self.xp.amax(evalf0, axis=tt)
#            f_min = self.xp.amin(evalf0, axis=tt)
#                        
#            y_rand = f_min + self.xp.random.random(N, dtype=self.evalType)*(f_max-f_min)
#            
#            ind_below = self.xp.where(y_rand < f(*integrationAndParamsVarSamplingGrids), axis=tt)
#
#            return f_max*(x1-x0)*(ind_below[0].shape[0])/N

    
    # implemented methods: rect, trap,
    def IntegralEval(self, lh, integral, paramsAndRanges, integrationVarsSamplingSchemes=None, method='trap'):        
        paramVariables = []
        paramNames, parametersLows, parametersHighs, parametersPoints, parametersSpacings = map(list, zip(*paramsAndRanges))        
        for paramName in paramNames:
            paramVariables.append(getSymbolByName(lh, paramName))
        integrationVariables, integrationVariablesLows, integrationVariablesHighs = map(list, zip(*integral.limits))          
        nVar = len(integrationVariables)        
        if integrationVarsSamplingSchemes is None:
            integrationVarsPoints =  [parametersPoints[0]] * nVar
            integrationVarsSpacings = ['linear'] * nVar
        else:
            if len(integrationVarsSamplingSchemes) == nVar:
                integrationVarsPoints, integrationVarsSpacings = map(list, zip(*integrationVarsSamplingSchemes))
            else:
                print(nVar)
                print(len(integrationVarsSamplingSchemes))
                print("Wrong sampling schemes for integration variables (check number of elements)")
                return
#        weights_S = sp.symbols('weights_S', real=True)
#        ff = weights_S * integral.function
        lambdaIntegrand = sp.lambdify((*integrationVariables, *paramVariables), integral.function, self.modules)

        integrationVariableSamplings = []
        smallIntegrationVariableSamplings = []
        parameterSamplings = []
        integrationVarsScheme = zip( integrationVariables, integrationVariablesLows, integrationVariablesHighs, integrationVarsPoints, integrationVarsSpacings)
        for integrationVariable, integrationVariableLow, integrationVariableHigh, integrationVarPoints, integrationVarSpacing in integrationVarsScheme:
            # print(integrationVariable, integrationVariableLow, integrationVariableHigh, integrationVarPoints, integrationVarSpacing)            
            if method=='rect':
                dx = (float(integrationVariableHigh)-float(integrationVariableLow))/float(integrationVarPoints-1)
                s = self.getSampling(float(integrationVariableLow)+dx/2, float(integrationVariableHigh)-dx/2, integrationVarPoints-1, integrationVarSpacing)
            elif method=='trap':
                s = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, integrationVarSpacing)
            elif method=='mc':
                s = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, 'random')
                s_small = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), max(100, integrationVarPoints/100), 'linear')
                smallIntegrationVariableSamplings.append(s_small)

            integrationVariableSamplings.append(s)

        # print(integrationVariableSamplings)
        for paramName, parameterLow, parameterHigh, parameterPoints, parameterSpacing in paramsAndRanges:    
            # print(paramName, parameterLow, parameterHigh, parameterPoints, parameterSpacing)
            parameterSamplings.append(self.getSampling(float(parameterLow), float(parameterHigh), parameterPoints, parameterSpacing))
        return self.outputData(parameterSamplings), self.outputData(self.intensityMap(self.parametricIntegralEvaluation(integrationVariableSamplings, parameterSamplings, lambdaIntegrand, method, smallIntegrationVariableSamplings)))
        

def polar_to_cart(polar_data, theta_step, range_step, x, y, order=3):
    from scipy.ndimage.interpolation import map_coordinates as mp
    X, Y = np.meshgrid(x, y)
    Tc = np.degrees(np.arctan2(Y, X)).ravel()
    Rc = (np.sqrt(X**2 + Y**2)).ravel()
    Tc[Tc < 0.0] = 360.0 + Tc[Tc < 0.0]
    Tc = Tc / theta_step
    Rc = Rc / range_step
    coords = np.vstack((Tc, Rc))
    polar_data = np.vstack((polar_data, polar_data[-1, :]))
    cart_data = mp(
        polar_data,
        coords,
        order=order,
        mode='constant',
        cval=np.nan)
    return(cart_data.reshape(len(y), len(x)).T)


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if a.dtype not in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array(a.shape)
    ndims = len(a.shape)
    if len(newdims) != ndims:
        print("[congrid] dimensions error. "
              "This routine currently only support "
              "rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray(newdims, dtype=float)
    dimlist = []

    if method == 'neighbour':
        for i in range(ndims):
            base = np.indices(newdims)[i]
            dimlist.append((old[i] - m1) / (newdims[i] - m1)
                           * (base + ofs) - ofs)
        cd = np.array(dimlist).round().astype(int)
        newa = a[list(cd)]
        return newa

    elif method in ['nearest', 'linear']:
        # calculate new dims
        for i in range(ndims):
            base = np.arange(newdims[i])
            dimlist.append((old[i] - m1) / (newdims[i] - m1)
                           * (base + ofs) - ofs)
        # specify old dims
        olddims = [np.arange(i, dtype=np.float) for i in list(a.shape)]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d(
            olddims[-1], a, kind=method, fill_value="extrapolate")
        newa = mint(dimlist[-1])

        trorder = [ndims - 1] + list(range(ndims - 1))
        for i in range(ndims - 2, -1, -1):
            newa = newa.transpose(trorder)

            mint = scipy.interpolate.interp1d(
                olddims[i], newa, kind=method, fill_value="extrapolate")
            newa = mint(dimlist[i])

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose(trorder)

        return newa
    elif method in ['spline']:
        oslices = [slice(0, j) for j in old]
        oldcoords = np.ogrid[oslices]
        nslices = [slice(0, j) for j in list(newdims)]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        # make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n",
              "Currently only \'neighbour\', \'nearest\',\'linear\',",
              "and \'spline\' are supported.")
        return None

from lambdifyDictionaries import *
from sympyHelpers import *

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
        self.modules = backendLibs[self.xp]        
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
            return l + self.xp.random.random(npoints)*(h-l)
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

    
    def parametricEvaluation(self, paramsSamplings, evalFunction):
        paramsSamplingGrids = self.xp.meshgrid( *paramsSamplings, sparse=True, copy=False)
        nParams = len(paramsSamplings)
        if self.xp == cp:
            @cp.fuse(kernel_name='evaluableFunction')
            def evaluableFunction(*paramsSamplingGrids):
                return evalFunction(*paramsSamplingGrids)
        else:
            evalFunctionV = np.vectorize(evalFunction)
            def evaluableFunction(*paramsSamplingGrids):
                return evaluFunctionV(*paramsSamplingGrids)
        
#        print(inspect.getmembers(evalFunction))
#        print(inspect.getsource(evalFunction))
#        print(inspect.getsource( inspect.unwrap(evaluableFunction)) )
        
        return evaluableFunction(*paramsSamplingGrids)


    def functionEval(self, functionToEval, paramsAndRanges):        
        paramVariables = []
        paramNames, parametersLows, parametersHighs, parametersPoints, parametersSpacings = map(list, zip(*paramsAndRanges))        
        for paramName in paramNames:
            paramVariables.append(getSymbolByName(functionToEval, paramName))
        nVar = len(paramVariables)
        lambdaIntegrand = sp.lambdify([*paramVariables], functionToEval, self.modules)
        parameterSamplings = []
        for paramName, parameterLow, parameterHigh, parameterPoints, parameterSpacing in paramsAndRanges:
            if parameterSpacing!='provided':
                parameterSamplings.append(self.getSampling(float(parameterLow), float(parameterHigh), parameterPoints, parameterSpacing))
            else:
                parameterSamplings.append(self.xp.asarray(parameterLow, dtype=self.evalType))
        return self.outputData(parameterSamplings), self.outputData(self.parametricEvaluation(parameterSamplings, lambdaIntegrand))

    
    def parametricIntegralEvaluation(self, integrationVarsSamplings, paramsSamplings, integrandFunction, method='rect', smallSamplings=[]):
        integrationAndParamsVarSamplingGrids = self.xp.meshgrid( *[*integrationVarsSamplings, *paramsSamplings], sparse=True, copy=False)
        weightsArrays = [self.getWeightsArray(aa.shape[0]) for aa in integrationVarsSamplings]
        # for Trapezoidal method only:
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

        # for Rectangles method
        def genericSum(xx):
            return self.xp.sum( xx, axis=tt )
            
        # for Trapezoidal method
        def genericReduction(xx):
            return self.postMap( self.xp.sum( xx * weightsGrid, axis=tt ) )

        N = 1
        for ssx in integrationVarsSamplings:
            N *= ssx.shape[0]
        M = paramsSamplings[0].shape[0]

        def mcReduction(y_comp):
            # integrationAndParamsVarSamplingGridsSmall = self.xp.meshgrid( *[*smallSamplings, *paramsSamplings], sparse=True, copy=False)
            # # evalf0 = integratedFunction(*integrationAndParamsVarSamplingGridsSmall)
            tt = (0,1)
            y_rand = self.xp.random.random( (integrationVarsSamplings[0].shape[0], integrationVarsSamplings[1].shape[0], M))
            f_range_real = self.xp.max(self.xp.real(y_comp), axis=tt)-self.xp.min(self.xp.real(y_comp), axis=tt)
            f_range_imag = self.xp.max(self.xp.imag(y_comp), axis=tt)-self.xp.min(self.xp.imag(y_comp), axis=tt)
            freqs_real = (self.xp.sum(self.xp.where(  y_rand*f_range_real < self.xp.real(y_comp) - self.xp.min(self.xp.real(y_comp), axis=tt), 1, 0 ), axis=tt) * f_range_real/float(N) + self.xp.min(self.xp.real(y_comp), axis=tt)) 
            freqs_imag = (self.xp.sum(self.xp.where(  y_rand*f_range_imag < self.xp.imag(y_comp) - self.xp.min(self.xp.imag(y_comp), axis=tt), 1, 0 ), axis=tt) * f_range_imag/float(N) + self.xp.min(self.xp.imag(y_comp), axis=tt))
            return freqs_real + 1j * freqs_imag

        
        if self.xp == cp:
            if method=='rect' or method=='raw':
                @cp.fuse(kernel_name='integratedFunctionRect')
                def integratedFunction(*integrationAndParamsVarSamplingGrids, reduce=genericSum, post_map=self.postMap):                               
                    return post_map(reduce( integrandFunction(*integrationAndParamsVarSamplingGrids))) 
            elif method=='trap':
                @cp.fuse(kernel_name='integratedFunctionTrap')
                def integratedFunction(*integrationAndParamsVarSamplingGrids):                               
                    return integrandFunction(*integrationAndParamsVarSamplingGrids)
            else:
                @cp.fuse(kernel_name='integratedFunctionMC')
                def integratedFunction(*integrationAndParamsVarSamplingGrids):
                    return integrandFunction(*integrationAndParamsVarSamplingGrids)
#                @cp.fuse(kernel_name='integratedFunctionMC')
#                def integratedFunction(*integrationAndParamsVarSamplingGrids, reduce=mcReduction, post_map=self.postMap):
#                    return integrandFunction(*integrationAndParamsVarSamplingGrids)
                
        else:
            integrandFunctionV = np.vectorize(integrandFunction)
            if method=='rect' or method=='raw':
                def integratedFunction(*integrationAndParamsVarSamplingGrids, post_map=self.postMap):
#                   return post_map(genericSum(np.nan_to_num(integrandFunction(*integrationAndParamsVarSamplingGrids))))
                    return post_map(genericSum(integrandFunctionV(*integrationAndParamsVarSamplingGrids)))
            else:
                def integratedFunction(*integrationAndParamsVarSamplingGrids):
                    return integrandFunctionV(*integrationAndParamsVarSamplingGrids)
            
        scaleFactor = 1.0
        for ii in range(nVars):
            if (len(integrationVarsSamplings[ii])>1):
#                scaleFactor *= (integrationVarsSamplings[ii][1] - integrationVarsSamplings[ii][0])
                scaleFactor *= (integrationVarsSamplings[ii][1] - integrationVarsSamplings[ii][0])
        if method=='trap':
            return scaleFactor * genericReduction(integratedFunction(*integrationAndParamsVarSamplingGrids))
        elif method=='rect' or method=='raw':
            return scaleFactor * integratedFunction(*integrationAndParamsVarSamplingGrids)
        elif method=='mc':        
            scaleFactor = 1.0
            tt = (0,1)
            for ii in range(nVars):
                if (len(smallSamplings[ii])>1):
                    scaleFactor *= (smallSamplings[ii][-1] - smallSamplings[ii][0])
            return self.xp.absolute(mcReduction(integratedFunction(*integrationAndParamsVarSamplingGrids)) * scaleFactor)
#            return integratedFunction(*integrationAndParamsVarSamplingGrids) * scaleFactor


    # implemented methods: rect, trap, mc
    def IntegralEval(self, lh, integral, paramsAndRanges, integrationVarsSamplingSchemes=None, method='rect'):        
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
            if method=='rect':
                dx = (float(integrationVariableHigh)-float(integrationVariableLow))/float(integrationVarPoints-1)
                s = self.getSampling(float(integrationVariableLow)+dx/2, float(integrationVariableHigh)-dx/2, integrationVarPoints-1, integrationVarSpacing)
            elif method=='raw':
                s = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, integrationVarSpacing)
            elif method=='trap':
                s = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, integrationVarSpacing)
            elif method=='mc':
                s = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, 'random')
                s_small = self.getSampling(float(integrationVariableLow), float(integrationVariableHigh), integrationVarPoints, 'linear')
                smallIntegrationVariableSamplings.append(s_small)
            else:
                print('Integration Method not valid or missing')
                return

            integrationVariableSamplings.append(s)
        for paramName, parameterLow, parameterHigh, parameterPoints, parameterSpacing in paramsAndRanges:    
            # print(paramName, parameterLow, parameterHigh, parameterPoints, parameterSpacing)
            if parameterSpacing!='provided':
                parameterSamplings.append(self.getSampling(float(parameterLow), float(parameterHigh), parameterPoints, parameterSpacing))
            else:
                parameterSamplings.append(self.xp.asarray(parameterLow, dtype=self.evalType))

        return self.outputData(parameterSamplings), self.outputData(self.intensityMap(self.parametricIntegralEvaluation(integrationVariableSamplings, parameterSamplings, lambdaIntegrand, method, smallIntegrationVariableSamplings)))

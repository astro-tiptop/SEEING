import matplotlib.pyplot as plt
import numpy as np
import scipy
from seeing.sympyHelpers import *
from sympy.parsing.sympy_parser import parse_expr

# in a formulary you annot have the same string as representation
# for two different symbols, representing two different entities
# if you put multiple formulas inside a dictionary, each defining 
# and using a symbol represented by an "f", in fact you are have multiple,
# unrelated symbols with the same representation. So a reworking of the formulas
# is needed to unify the symbols variable in order to have a common "f" symbols
# shared by the formulas. This allows a reasonable composition of the expressions.


class Formulary(object):
    # each element of the formulas list is a tuple (_lhs, _rhs, _eq)
    def __init__(self, name='', names=[], formulas=[], symbol_map={}):
        self.name = name
        self.symbols = set({})
        self.functions = set({})
        self.formulas = dict(zip(names, formulas))
        self.symbol_map = symbol_map
        for ff in formulas:
            if type(ff) is tuple:
                for ffs in ff:
                    self.symbols.update(ffs.free_symbols)
                    if type(ffs)==sp.Eq:
                        self.functions.update(ffs.lhs.atoms(sp.Function))
                    else:
                        self.functions.update(ffs.atoms(sp.Function))                    
            else:
                self.symbols.update(ff.free_symbols)
                if type(ff)==sp.Eq:
                    self.functions.update(ff.lhs.atoms(sp.Function))
                else:
                    self.functions.update(ff.atoms(sp.Function))                    


    def writeToFile(self, filename):
        f = open(filename, "w")
        for name, expr in self.formulas.items():
            f.write(name+"\n")
            f.write(str(expr)+"\n")
        f.close()
        

    def loadFromFile(filename):
        names = []
        formulas = []
        f = open(filename, "r")
        for aname in f:
            names.append(aname.replace('\n', ''))
            ss = f.readline().replace('\n', '')
            ss = ss.replace('lambda', 'lambda_')
            ss = ss.replace('^', '__')
            aformula = parse_expr(ss, evaluate=False)
            aformula = subsParamsByName(aformula, {'lambda_': sp.symbols('lambda')})
            
            aformula = subsSpecialParams(aformula, '__', '^')

            formulas.append(aformula)
        f.close()
        return Formulary('A', names, formulas)
        
        
    def addFormula(self, name, formula):
        self.formulas[name] = formula

    def displaySymbols(self):
        display(sp.Matrix([[*self.symbols]]))

    def displayFunctions(self):
        display(sp.Matrix([[*self.functions]]))
        
    def displayAll(self, showNames=True):
        for name, f in self.formulas.items():
            if showNames:
                print(name)
            display(f)

    def display(self, aa, showNames=True):
        for name, f in self.formulas.items():
            if name in aa:
                if showNames:
                    print(name)
                display(f)
            
    def getFormula(self, name):
        return self.formulas[name]

    def __getitem__(self, key):
        return self.formulas[key]
    
    def getFormulaRhs(self, name):
        return self.formulas[name].rhs

    def getFormulaLhs(self, name):
        return self.formulas[name].lhs

    def getRestrictedLambda(
            self,
            name,
            subsDict,
            independentVarNames,
            modules=cpulib):

        _expr = self.getFormula(name)
        if type(_expr)==sp.Eq:
            _expr = _expr.rhs
        return getRestrictedLambda(
            _expr, subsDict, independentVarNames, modules)

    def evaluateFormula(
            self,
            name,
            subsDict,
            samples,
            mCalc):
        _expr = self.getFormula(name)
        if type(_expr)==sp.Eq:
            _expr = _expr.rhs
        _expr = subsParamsByName(_expr, subsDict)
        return mCalc.functionEval(_expr, samples)

    def plotFormula(
            self,
            name,
            subsDict,
            samples,
            mCalc,
            log=False):
        _params, _data = self.evaluateFormula(
            name, subsDict, samples, mCalc)
        nvars = len(samples)
        if nvars == 1:
            fig = plt.figure(figsize=(7, 7))
            axs = fig.add_subplot(111)
            axs.axis('auto')
            if log:
                plt.xscale('log')
                plt.yscale('log')
            axs.plot(_params[0], np.real(_data))
        elif nvars == 2:
            fig, ax = plt.subplots(figsize=(7,7))
            if log:
                ax.imshow(np.log(np.abs(_data)), cmap='hot')
            else:
                ax.imshow(np.real(_data), cmap='hot')
        else:
            pass
        plt.show()

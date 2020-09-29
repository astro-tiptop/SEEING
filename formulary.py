import matplotlib.pyplot as plt
import numpy as np
import scipy
from sympyHelpers import *

# in a formulary you annot have the same string as representation
# for two different symbols, representing two different entities
# if you put multiple formulas inside a dictionary, each defining 
# and using a symbol represented by an "f", in fact you are have multiple,
# unrelated symbols with the same representation. So a reworking of the formulas
# is needed to unify the symbols variable in order to have a common "f" symbols
# shared by the formulas. This allows a reasonable composition of the expressions.


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
        _expr = self.getFormula(name).rhs
        return getRestrictedLambdaBasic(
            _expr, subsDict, independentVarNames, modules)

    def evaluateFormula(
            self,
            name,
            subsDict,
            independentVarNames,
            samples,
            modules=cpulib):
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

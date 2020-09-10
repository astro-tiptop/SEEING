import sympy as sp
import inspect
from functools import reduce
from SEEING.lambdifyDictionaries import *


def getSymbolByName(expr, sname):
    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
    if sname in free_syms.keys():
        return free_syms[sname]
    else:
        return None
    

def subsParamsByName(expr1, subsDict):
    expr = expr1
    for sname, sval in subsDict.items():
        ss = getSymbolByName(expr, sname)
        if ss:
            expr = expr.subs(ss, sval)
    return expr


def getRestrictedLambdaBasic(
        expr,
        subsDict,
        independentVarNames,
        lmodules=cpulib):
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
    return sp.lambdify(tt, expr, lmodules)


def evaluateLambda(p_f, independentVarNames, samples, xp=np):
    nvars = len(independentVarNames)
    if nvars==0:
        return None
    elif nvars==1 or samples[0].ndim != 1:
        return (*samples, p_f(*samples))
    elif nvars>1 and samples[0].ndim == 1:
        varsSamplingGrids = xp.meshgrid( *samples, sparse=False, copy=True)
        return (*varsSamplingGrids, p_f(*varsSamplingGrids))    
    else:
        return None


def evaluateFormula(
        expr,
        subsDict,
        independentVarNames,
        samples,
        modules=cpulib):
    p_f = getRestrictedLambdaBasic(expr, subsDict, independentVarNames, modules)
    return evaluateLambda(p_f, independentVarNames, samples)


from sympy.parsing import latex

def latexToSympy(latexText):
    return latex.parse_latex(latexText)


def sympyToLatex(expr):
    return sp.latex(expr, mode='plain')
    

def sympyToString(expr):
    return sp.srepr(expr)
    
    
    
    
# Currently not used
#def replaceParamsByName(expr, subsDict):    
#    syms_and_funs = set(expr.free_symbols) | set([i.func for i in expr.atoms(sp.Function) if isinstance(i, sp.function.AppliedUndef)])
#    free_syms = {symbol.name: symbol for symbol in syms_and_funs}
#    for sname, sval in subsDict.items():
#        print(sname, sval)
#        expr = expr.replace({free_syms[sname], sval})
#    return expr


#ssr = sp.srepr(expr_Phi)
#print(ssr)
#from sympy.parsing.sympy_parser import parse_expr
#
#ssrE = parse_expr(ssr)
#display(ssrE)

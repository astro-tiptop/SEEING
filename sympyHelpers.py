import sympy as sp
import inspect
from functools import reduce
from lambdifyDictionaries import *


def getSymbolByName(expr, sname):
    """Gets a symbol from a sympy expressions, given the symbols name

    Parameters
    ----------
    expr : sympy expression
        A sympy expression 
    sname : string
        A string containing the symbol name

    Returns
    -------
    Symbol
        The symbol having the required name. None is returned if the symbol is not found
    """

    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}
    if sname in free_syms.keys():
        return free_syms[sname]
    else:
        return None
    

def subsParamsByName(expr, subsDict):
    """Generate a new sympy expression from a given one by substitution of a set of its symbols. 
    Substitutions are specified as a dictionay of string:value pairs (value can be any sympy expression 
    or numerical value).

    Parameters
    ----------
    expr : sympy expression
        A sympy expression 
    subsDict : dict
        A dictionary containing the "symbol name":value substitution pairs

    Returns
    -------
    expr
        The resulting sympy expression
    """    

    for sname, sval in subsDict.items():
        ss = getSymbolByName(expr, sname)
        if ss:
            expr = expr.subs(ss, sval)
    return expr


def lambdifyByName(expr, symbolsNames, modules_):
    """Wraps around sympy.lambdify to generate a lambda function from a given sympy expression. 
    The set of symbols to be passed to sympy.lambdify is specified as a list of strings. 

    Parameters
    ----------
    expr : sympy expression
        A sympy expression 
    symbolsNames : string
        A list of strings containing the symbols name
    modules_ : see sympy lambdify modules 

    Returns
    -------
    lambda
        The resulting lambda 
    """    

    symbolsList = []
    for sName in symbolsNames:
        symbolsList.append(getSymbolByName(expr, sName))
    return sp.lambdify( symbolsList, expr, modules=modules_ )


def getRestrictedLambda(
        expr,
        subsDict,
        independentVarNames,
        modules=cpulib):
    """Wraps around sympy.lambdify to generate a lambda function from a given sympy expression. 
    The set of symbols to be passed to sympy.lambdify is specified as a list of strings. 
    If one or more symbols are not found in the expressions, dummy symbols are introduced.

    Parameters
    ----------
    expr : sympy expression
        A sympy expression 
    symbolsNames : string
        A list of strings containing the symbols names
    modules_ : see sympy lambdify modules 

    Returns
    -------
    lambda
        The resulting lambda 
    """    

    
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
    return sp.lambdify(tt, expr, modules)


def evaluateLambda(p_f, independentVarNames, samples, xp=np):
    """Evaluates a given a ufunc (numpy or cupy)

    Parameters
    ----------
    p_f : ufunc
        A numpy or cupy ufunc 
    independentVarNames : string list 
        A list of strings containing the names of the symbols to be used as evaluations variables
    samples : a list of numpy or cupy arrays 
        A samples vector for each variable, to be used generate the grid to evaluate the ufunc
    xp : Array library, either np (numpy) or cp (cupy)

    Returns
    -------
    array
        The result of the function evaluation 
    """

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
    """Evaluates a given a sympy expression

    Parameters
    ----------
    expr : sympy expressions
        The expressions to be evaluated
    subsDict : dict
        A dictionary containing the "symbol name":value substitution pairs
    independentVarNames : string list 
        A list of strings containing the names of the symbols to be used as evaluations variables
    samples : a list of numpy or cupy arrays 
        A samples vector for each variable, to be used generate the grid to evaluate the expressions
    modules : either cpulib (numpy) or gpulib (cupy)

    Returns
    -------
    array
        The result of the expressions evaluation 
    """

    p_f = getRestrictedLambdac(expr, subsDict, independentVarNames, modules)
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


def subsSpecialParams(expr, s1, s2):
    
    free_syms = {symbol.name: symbol for symbol in expr.free_symbols}

    for sname in free_syms:
        ss = getSymbolByName(expr, sname)
        expr = expr.subs(ss, sp.symbols(sname.replace(s1, s2)))

    return expr

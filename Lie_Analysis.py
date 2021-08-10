import copy

import pandas as pd

from utils import DEint as de

from sympy import *

import regex as re

# RMHD_CHM
# var_list = ['x', 'y', 'z', 't', 'tau', 'phi', 'psi', 'chi']
# var_dict = {'xse1':'xi^(1)', 'xse2':'xi^(2)', 'xse3':'xi^(3)', 'xse4':'xi^(4)',
#              'xse5':'xi^(5)', 'eta1':'eta^(1)', 'eta2':'eta^(2)', 'eta3':'eta^(3)'}


def Symmetry_Analysis(cvs_path, N_indep, N_dep,
                      var_list, simplify=True, latex=False):
    """Gives a symbolic representation of all equations
       from a cvs file.

       Args:
       cvs_path (str): path where the cvs is.
       N_indep (int): number of independent variables.
       N_dep (int): number of dependant variables.
       var_list (list): list of all variables labels
                        in string format.
       simplify (boolean): If True, simplifies the set
                           of determining equations.
    """
    var_dict = {}
    for i in range(0, N_indep):
        var_dict['xse' + str(i+1)] = 'xi^' + '(' + var_list[i] + ')'
    for i in range(0, N_dep):
        k = i + N_indep
        if len(var_list[k]) > 1:
            var_dict['eta' + str(i+1)] = 'eta^' + \
                '(' + "\\" + var_list[k] + ')'
        else:
            var_dict['eta' + str(i+1)] = 'eta^' + "(" + var_list[k] + ')'
    latex_dict = {}
    for i in range(0, N_indep):
        if len(var_list[i]) > 1:
            latex_dict['xse' + str(i+1)] = 'xi^' + \
                '{' + "\\" + var_list[i] + '}'
        else:
            latex_dict['xse' + str(i+1)] = 'xi^' + "{" + var_list[i] + '}'
    for i in range(0, N_dep):
        k = i + N_indep
        if len(var_list[k]) > 1:
            latex_dict['eta' + str(i+1)] = 'eta^' + \
                '{' + "\\" + var_list[k] + '}'
        else:
            latex_dict['eta' + str(i+1)] = 'eta^' + "{" + var_list[k] + '}'
    det_eqns, constants = cvs_to_list(cvs_path, var_list)
    if simplify:
        det_eqns = simplify_redundant_eqn(det_eqns)
        det_eqns = simplify_redundant_eqn_second_phase(det_eqns)
    if latex:
        latex_code = latex_det_eqn(det_eqns,
                        latex_dict, var_list, constants)
        latex_code = latex_code.replace("+-", "-")
        return print(latex_code)
    return sym_det_eqn(det_eqns, var_dict, var_list, constants)


def cvs_to_list(cvs_path, var_list):
    """Converts a cvs containing all the determining equiation
       into a list of lists. Each internal list will have a
       dictionary containing all the relevant information
       of each term.

       cvs_path (str): path where the cvs is.
    """
    raw_eqs = pd.read_csv(cvs_path,  header=None)
    pd.set_option('display.max_colwidth', None)
    raw_eqs.columns = ['eqn']
    raw_eqs.eqn = raw_eqs.eqn.apply(lambda eq: eq.strip(' == 0'))
    constants = de.find_constants(raw_eqs.eqn, var_list)
    det_eqn = {}
    for i, row in raw_eqs.iterrows():
        det_eqn[i] = de.eqn_process(row.eqn, constants, var_list)
    return det_eqn, constants


def simplify_redundant_eqn(det_eqns):
    """given a dict of equations reduces the set 
       by eliminating redundant eqns. It is just a test
       to try the logic far from being ready yet

       Args:
       det_eqns (dict): 
    """
    simplify = True
    N = len(det_eqns)
    exit = 0
    zero_terms = {}
    det_eqns_aux = copy.deepcopy(det_eqns)
    while simplify:
        for idx, eqn in det_eqns.items():
            for i, zero in zero_terms.items():
                for term in eqn:
                    if de.is_zero(zero, term) and term in det_eqns_aux[idx]:
                        det_eqns_aux[idx][det_eqns_aux[idx].index(term)] = 0
                        exit = 0
                det_eqns_aux[idx] = list(
                    filter(lambda num: num != 0, det_eqns_aux[idx]))
                exit += 1
            if len(det_eqns_aux[idx]) == 1:
                zero_terms[idx] = det_eqns_aux[idx][0]
            if exit > N:
                simplify = False
        det_eqns = copy.deepcopy(det_eqns_aux)
    simplify_det_eqns = {k: v for k, v in det_eqns.items() if v}
    for idx in zero_terms:
        simplify_det_eqns[idx] = [zero_terms[idx]]
    return simplify_det_eqns


def simplify_redundant_eqn_second_phase(det_eqn):
    for idx, eqn in det_eqn.items():
        det_eqn[idx] = de.drop_constants(eqn)
    return det_eqn


def latex_eqn_code(eqn, var_dict, var_list, constants):
    """This code generates the latex string format of 
    the determinant equations

    Parameters
    ----------
    eqn : [dictionary]
        dictionary containing all the determinant eqautions
    var_dict : [dictionary]
        dictionary mapping the name of the variable to it's 
        latex code 
    var_list : [list]
        list with all independant and dependant variables in
        string format
    constants : [list]
        list containing all constants

    Returns
    -------
    [string]
        returns a single string with  all the equations 
        typed in latex format
    """
    sym_cte_list = []
    for idx in range(len(var_list)):
        sym_cte_list.append(var_list[idx] + str(idx))
    for idx in range(len(constants) - len(var_list)):
        sym_cte_list.append('alpha_' + str(idx))
    A = ''
    for term in eqn:
        one_term = False
        if len(eqn) == 1:
            one_term = True
        A = A + de.dict_to_latex(term, var_dict, var_list,
                                 sym_cte_list, one_term) + '+'
    A = A[:-1]
    A = A + '=0'
    return A


def latex_det_eqn(det_eqn, var_dict, var_list, constants):
    """Gives the latex code version of remaining
       determining equation.

       Args:
       det_eqn (dict): dictionary will all the
                       determining equations. 
    """
    latex_code = ''
    for eqn in det_eqn.values():
        latex_code = latex_code + latex_eqn_code(eqn,
                                                 var_dict, var_list, constants) + "\\" + "\\" + "\n"
    return latex_code


def get_symbolic_terms(eqn, var_dict, var_list, constants):
    """given a list of dictionaries with the
       information of each term, retuns the 
       symbolic equivalent.

       Args:
       eqn (list): list of dictionaries 
    """
    sym_cte_list = []
    for idx in range(len(var_list)):
        sym_cte_list.append(symbols(var_list[idx] + str(idx)))
    for idx in range(len(constants) - len(var_list)):
        sym_cte_list.append(symbols('alpha_' + str(idx)))
    A = 0
    for term in eqn:
        one_term = False
        if len(eqn) == 1:
            one_term = True
        A += de.dict_to_symb(term, var_dict, var_list,
                             sym_cte_list, one_term)
    return A


def sym_det_eqn(det_eqn, var_dict, var_list, constants):
    """Gives the symbolic version of remaining
       determining equation.

       Args:
       det_eqn (dict): dictionary will all the
                       determining equations. 
    """
    M = Matrix([[]])
    i = 1
    for eqn in det_eqn.values():
        M = M.row_insert(i-1, Matrix([[i,
                                       Eq(get_symbolic_terms(
                                           eqn, var_dict, var_list, constants), 0)
                                       ]]))
        i += 1
    return M

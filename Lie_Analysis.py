import copy

import pandas as pd

from utils import DEint as de

from sympy import *

# Global variables
x, y, z, t, tau = symbols('x y z t tau')
phi, psi, chi = symbols('phi psi chi')

xi_1 = Function('xi_1')(x, y, z, t, tau, phi, psi, chi)
xi_2 = Function('xi_2')(x, y, z, t, tau, phi, psi, chi)
xi_3 = Function('xi_3')(x, y, z, t, tau, phi, psi, chi)
xi_4 = Function('xi_4')(x, y, z, t, tau, phi, psi, chi)
xi_5 = Function('xi_5')(x, y, z, t, tau, phi, psi, chi)

eta_1 = Function('eta_1')(x, y, z, t, tau, phi, psi, chi)
eta_2 = Function('eta_2')(x, y, z, t, tau, phi, psi, chi)
eta_3 = Function('eta_3')(x, y, z, t, tau, phi, psi, chi)

var_dict = {'xse1':xi_1, 'xse2':xi_2, 'xse3':xi_3, 'xse4':xi_4, 'xse5':xi_5, 
           'eta1':eta_1, 'eta2':eta_2, 'eta3':eta_3}
var_list = [x, y, z, t, tau, phi, psi, chi] 

def cvs_to_list():
    """Converts a cvs containing all the determining equiation
       into a list of lists. Each internal list will have a
       dictionary containing all the relevant information
       of each term.

    """
    raw_eqs = pd.read_csv('RMHD_CHM_DE.csv', header=None)
    pd.set_option('display.max_colwidth', None)
    raw_eqs.columns = ['eqn']
    raw_eqs.eqn = raw_eqs.eqn.apply(lambda eq: eq.strip(' == 0'))
    global constants
    constants = de.find_constants(raw_eqs.eqn)
    det_eqn = {}
    for i, row in raw_eqs.iterrows():
        det_eqn[i] = de.eqn_process(row.eqn, constants)
    return det_eqn

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
    det_eqns_aux =copy.deepcopy(det_eqns)
    while simplify:
        for idx, eqn in det_eqns.items():
            for i, zero in zero_terms.items():
                for term in eqn:
                    if de.is_zero(zero, term) and term in det_eqns_aux[idx]:
                        det_eqns_aux[idx][det_eqns_aux[idx].index(term)] = 0                       
                        exit = 0
                det_eqns_aux[idx] = list(filter(lambda num: num != 0, det_eqns_aux[idx]))
                exit += 1
            if len(det_eqns_aux[idx]) == 1:
                zero_terms[idx] = det_eqns_aux[idx][0]
            if exit > N:
                simplify = False
        det_eqns =copy.deepcopy(det_eqns_aux)
    simplify_det_eqns = {k:v for k,v in det_eqns.items() if v}
    for idx in zero_terms:
        simplify_det_eqns[idx] =[zero_terms[idx]]
    return simplify_det_eqns

def get_symbolic_terms(eqn):
    """given a list of dictionaries with the
       information of each term, retuns the 
       symbolic equivalent.

       Args:
       eqn (list): list of dictionaries 
    """
    sym_cte_list = []
    for idx in range(len(constants)):
        sym_cte_list.append(symbols('alpha_' + str(idx)))
    A = 0 
    for term in eqn:
        one_term = False
        if len(eqn) == 1:
            one_term = True
        A += de.dict_to_symb(term, var_dict, 
                            var_list, sym_cte_list, one_term)
    return A

def sym_det_eqn(det_eqn):
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
            get_symbolic_terms(eqn)
        ]]))
        i += 1
    return M
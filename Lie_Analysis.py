import copy

import pandas as pd

from utils import DEint as de


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
        print(exit)
        det_eqns =copy.deepcopy(det_eqns_aux)
    simplify_det_eqns = {k:v for k,v in det_eqns.items() if v}
    return simplify_det_eqns

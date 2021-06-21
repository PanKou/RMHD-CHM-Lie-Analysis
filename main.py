import Lie_Analysis as LA

det_eqns = LA.cvs_to_list()
print(LA.simplify_redundant_eqn(det_eqns))
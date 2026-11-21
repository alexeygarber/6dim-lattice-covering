import cvxpy as cp
import numpy as np
import warnings

import simplex
from fractions import Fraction as F

# Compute the matrix B_L using the input simplex V (rows are spanning vectors) and variable PQF Q

def B_L(V,Q):
    G = (V @ Q) @ V.T
    d = cp.reshape(cp.diag(G), (6, 1), order="F")
    B = cp.bmat([[np.array([[4.0]]), d.T], [d, G]])
    return 0.5*(B+B.T)


# simplices has to be a list of Simplex
# Output uses float and may be optimal_inaccurate
# We need rational confirmation to test the bound. This is done using rational_psd function

def maxdet(simplices):
    if len(simplices)==0:
        raise ValueError("List of simplices has to be nonempty")
    
    for i in simplices:
        if not isinstance(i, simplex.Simplex):
            raise ValueError("Only Simplex allowed in maxdet")
   
    Q = cp.Variable((6, 6), symmetric=True)
    constraints = [B_L(np.array(i.V, dtype=float),Q) >> 0 for i in simplices]

    obj = cp.Maximize(cp.log_det(Q))
    prob = cp.Problem(obj, constraints)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Solution may be inaccurate")
        try:
            prob.solve(solver=cp.CLARABEL)
        except cp.error.SolverError:
            return("solver_error",None,None,None)
    if prob.status == cp.OPTIMAL or prob.status == cp.OPTIMAL_INACCURATE:
        return (prob.status, float(np.exp(prob.value)), [i.dual_value for i in constraints],Q.value)
    else:
        return (prob.status,None,None,None)


# Approximate float matrix M with a positive semidefinite rational with denominator den
# This is used to get rational certificate from float dual values for a convex optimization problem

def rational_psd(M,den):
    if int(den)!=den:
        raise ValueError("Denominator should be an integer")
    den = int(den)
    if den <= 0:
        raise ValueError("Denominator should be positive")
    M = np.array(M, dtype=float)
    
    if M.shape != (7, 7):
        raise ValueError("The dual value should be 7x7 matrix")

    M = 0.5 * (M+M.T)
    lam, U = np.linalg.eigh(M)

    tmp = [[F(0,1)] * 7 for _ in range(7)]


    for k in range(7):
        mu = F(int(round(float(lam[k]) * den)), den)
        if mu<=0:
            continue
        W = [F(int(round(float(U[i][k]) * den)), den) for i in range(7)]
        for j in range(7):
            for l in range(7):
                tmp[j][l] += mu*W[j]*W[l]
    return tuple([tuple(i) for i in tmp])


    

# Testing

if __name__ == "__main__":
    
    L1=simplex.Simplex([[0,0,0,0,0,0],
                [1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]])

    L2=simplex.Simplex([[2,2,1,-1,-1,-1],
                [1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]])


    print(maxdet([L1,L2]))

    

    
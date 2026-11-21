import numpy as np
from pathlib import Path
import sys
import os
import json

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))

from edges import edge_candidates
from cutoff import EDGE_CUTOFF

def apexes_sigma0(edges_cand):
    """
    Function finds potential apexes for neighbors of the standard simplex across the facet x_1+...+x_6=1
    edges_cand is the list of admissible edges
    """

    for edge in edges_cand:
        if len(edge)!=6:
            raise ValueError("Each edge candidate must have 6 coordinates")
    
    apexes = []

    for edge in edges_cand:
        for i in [-1,1]:
            apex = np.array([1,0,0,0,0,0])+i*np.array(edge)
            if sum(apex)==2:
                tmp = True
                for j in range(6):
                    ej=np.array([0,0,0,0,0,0])
                    ej[j]=1
                    if (tuple(apex-ej) not in edges_cand) and (tuple(ej-apex) not in edges_cand):
                        tmp=False
                        break
                if tmp:
                    apexes.append(tuple([int(j) for j in apex]))
    return frozenset(apexes)

#
# For simplex of normalized volume 1, the smallest eigenvalue of the inertia matrix M_L is lambda > 0.128040
# 4/lambda < 32 so 31 is the cutoff we need for potential edges
# If x = sum x_k e_k , then barycentric coordinate of x are (-sum x_k, x_1, ..., x_6) 
# sum x_k^2 + (sum x_k)^2 <= 4/lambda
#

cand = apexes_sigma0(edge_candidates(EDGE_CUTOFF))

types = {}

for i in cand:
    tmp = tuple(sorted(i, reverse=True))
    if tmp in types:
        types[tmp] += 1
    else:
        types[tmp] = 1


script_dir = Path(__file__).resolve().parent
file_path = script_dir / "step2_apexes.json"

with open(file_path, "w") as f:
    tmp = sorted(types)
    f.write("[\n" + ",\n".join(json.dumps(t) for t in tmp) + "\n]\n")

with open(file_path, "r") as f:
    data = json.load(f)
    if ({tuple(t) for t in data}!=set(types)):
        raise ValueError("Writing to file was not successful")

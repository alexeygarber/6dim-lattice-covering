# Setting
# We start from the standard simplex sigma0
# Suppose its neighbor across the facet with vertices e1,...,e6 is of type A with apex [1,1,1,-1,0,0]
# Claim 1: Then the neighbor across the facet 0, e1, e2, e3, e5, e6 MUST have the same apex [1,1,1,-1,0,0] (and be of type A)

# The only assumption is that all neighbors are of types A, B, C


from pathlib import Path
import sys
import os
import json
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))

from simplex import Simplex
from relax import maxdet
from cutoff import EDGE_CUTOFF, D0
from edges import edge_candidates
from configuration import apex_coefficients, neighbor_candidates
from search import close_tree, process_child
from symmetry import canonical_form

total_start = datetime.now()

script_dir = Path(__file__).resolve().parent
with open(script_dir / "../02_step2_types/step2_ABC.json") as f:
    
    apexes = [[int(x) for x in a] for a in json.load(f)]

if apexes != [[1, 1, 0, 0, 0, 0], [1, 1, 1, 0, 0, -1], [1, 1, 1, 1, -1, -1]]:
    raise ValueError("Incorrect types read")

apexes = tuple(tuple(apexes[i]) for i in [1,0,2])

edges = edge_candidates(EDGE_CUTOFF)
bary = apex_coefficients(apexes)


# standard simplex
sigma0=Simplex([[0,0,0,0,0,0],
                [1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]])

# new apex
new_apex = [1,1,1,-1,0,0]

# A-type neighbor
sigma1=sigma0.replace([0,0,0,0,0,0],new_apex)

# We aim to prove that the target simplex is also present
target = sigma0.replace([0,0,0,1,0,0],new_apex)

det = maxdet([sigma0,sigma1])

if det[1] is None:
    raise ValueError("maxdet failed on type A neighbors")

if det[1]<float(D0):
    raise ValueError("That determinant should have been rejected by maxdet on two simplices which contradicts the type A structure")


# Creating list of candidates for the other facet
candidates = neighbor_candidates([sigma0, sigma1], sigma0.facet(4), [0,0,0,1,0,0], edges, bary)


# Verifying that the target simplex is in that list
if target.key in [i.key for i in candidates]:
    deviations = [i for i in candidates if i.key!=target.key]
else:
    raise ValueError("The target simplex is not in the list of candidates")

open_trees = []

closed = {
    "frame":0,
    "open":0,
    "certified":0
}


# Checking all candidates and keeping only those that are not refuted by frame and not certified
for i in deviations:
    tmp = process_child([sigma0,sigma1],i)
    closed[tmp[0]]+=1
    if tmp[0] == "open":
        tmp2=maxdet(tmp[1])[1]
        if tmp2 is None:
            raise ValueError("maxdet failed for a configuration")
        open_trees.append([tmp2,tmp[1],tmp[2]])

if len(deviations) != sum(closed[i] for i in closed):
    raise ValueError("Every configuration must be processed")


print("Total candidates:",len(candidates))
print("Total candidates other than target:",len(deviations))
print("Total closed by frame:",closed["frame"])
print("Total certified:",closed["certified"])
print("Total open remaining:",closed["open"],"\n")

# Creating dictionary of remaining configuration up to symmetry
tree_dict = {}

for i in open_trees:
    tmp = canonical_form(i[1])
    if tmp not in tree_dict:
        tree_dict[tmp] = i+[1]
    else:
        tree_dict[tmp][3]+= 1

cases = sorted(tuple(tree_dict[i] for i in tree_dict), key=lambda c: c[0])

print ("Total cases",len(cases))
for i in range(len(cases)):
    print("Case",i,"maxdet = ",cases[i][0],"vertex = ",cases[i][1][-1].pts[0],"Chosen from",cases[i][3],"options")

print("\n")




for i in cases:
    start = datetime.now()

    print("Working on vertex",i[1][-1].pts[0],".",i[3],"copies. maxdet =",i[0])
    tmp = close_tree([[i[1],i[2]]],edges,bary,100)
    
    end = datetime.now()
        
    print("Total visited for",i[1][-1].pts[0],":",tmp["total"])
    print("Total certified:",tmp["certified"])
    print("Total children:",tmp["children"])
    print("Total excluded by frame cutoff:",tmp["frame"])
    #print("Total excluded by library:",tmp["library"])
    print("Total duplicates:",tmp["duplicate"])
    #print("Total duplicate roots:",tmp["duplicate_roots"])
    print("Maximal configuration size:",tmp["max_size"])
    print("Total time spent: ", (end-start).total_seconds(),"\n")


print("Total time:", (datetime.now()-total_start).total_seconds(),"seconds")


    
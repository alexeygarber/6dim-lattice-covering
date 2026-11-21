# Setting
# We start from the standard simplex sigma0
# Suppose its neighbor sigma1 across the facet with vertices e1,...,e6 is of type A with apex [1,1,1,-1,0,0]
# Then Claims 1 and 2 give a collection of 9 simplices forming a free sum of three triangles
# Triangle 1: e1, e2, e3
# Triangle 2: 0, e4, [1,1,1,-1,0,0]
# Triangle 3: e5, e6, [1,1,1,0,-1,-1]

# The other three neighbors of sigma0 are of types B or C
# Let sigma be the neighbor across facet opposite to e1

# Claim: the free-sum cell of sigma does not depend on the choice of sigma and the other three vertices of this cell are
# (-1,0,0,0,1,1), (-1,0,0,1,0,0), (-1,-1,-1,1,1,1)

# The only assumption is that all neighbors are of types A, B, C


from pathlib import Path
import sys
import os
import json
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))

from simplex import Simplex, symmetry_key
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


# apexes are now sorted type A, type B, type C
apexes = tuple(tuple(apexes[i]) for i in [1,0,2])

edges = edge_candidates(EDGE_CUTOFF)
bary = apex_coefficients(apexes)
bary_BC = apex_coefficients([apexes[1],apexes[2]])

# standard simplex
sigma0=Simplex([[0,0,0,0,0,0],
                [1,0,0,0,0,0],
                [0,1,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,1,0,0],
                [0,0,0,0,1,0],
                [0,0,0,0,0,1]])

# new apex
type_A_apex = [1,1,1,-1,0,0]

# A-type neighbor
sigma1=sigma0.replace([0,0,0,0,0,0],type_A_apex)

# Another type A guaranteed by Claim 1
sigma2 = sigma0.replace([0,0,0,1,0,0],type_A_apex)

# This is the second type A neighbor guaranteed by Claim 2 and its companion from Claim 1
second_type_A_apex = [1,1,1,0,-1,-1]
sigma3 = sigma0.replace([0,0,0,0,1,0],second_type_A_apex)
sigma4 = sigma0.replace([0,0,0,0,0,1],second_type_A_apex)

# More type A neighbors of other simplices guaranteed by Claims 2 and 1
sigma5 = sigma1.replace([0,0,0,0,0,1],second_type_A_apex)
sigma6 = sigma1.replace([0,0,0,0,1,0],second_type_A_apex)
sigma7 = sigma2.replace([0,0,0,0,0,1],second_type_A_apex)
sigma8 = sigma2.replace([0,0,0,0,1,0],second_type_A_apex)

conf = [sigma0,sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8]

det = maxdet(conf)

if det[1] is None:
    raise ValueError("maxdet failed on a triangulation of the free sum")

if det[1]<float(D0):
    raise ValueError("That determinant should have been rejected by maxdet on nine simplices which contradicts existence of Q_6^c1")


# Creating list of candidates for the facet of sigma0 opposite to e1. Not using type A due to Claim 3
candidates = neighbor_candidates(conf, sigma0.facet(1), [1,0,0,0,0,0], edges, bary_BC)


open_trees = []

closed = {
    "frame":0,
    "open":0,
    "certified":0
}


# Checking all candidates and keeping only those that are not refuted by frame and not certified
for i in candidates:
    tmp = process_child(conf,i)
    closed[tmp[0]]+=1
    if tmp[0] == "open":
        tmp2=maxdet(tmp[1])[1]
        if tmp2 is None:
            raise ValueError("maxdet failed for a configuration")
        open_trees.append([tmp2,tmp[1],tmp[2]])

if len(candidates) != sum(closed[i] for i in closed):
    raise ValueError("Every configuration must be processed")

print("Total candidates:",len(candidates))
print("Total closed by frame:",closed["frame"])
print("Total certified:",closed["certified"])
print("Total open remaining:",closed["open"],"\n")

# Verifying that we have only three options remaining.
# All three options are simplices in triangulation of the free sum of three triangles
# 0, e4, (-1,0,0,0,1,1) 
# e2, e3, (-1,-1,-1,1,1,1)
# e5, e6, (-1,0,0,1,0,0) 
if len(open_trees)==3:
    print("The only potential neighbors of sigma0 across the facet opposite to (1,0,0,0,0,0) are")
    for i in open_trees:
        print("Simplex with vertices",i[1][-1].pts)
else:
    raise ValueError("We should have only three options remaining")


#raise ValueError("STOP")

# Creating dictionary of remaining configuration up to symmetries keeping sigma0
tree_dict = {}

for i in open_trees:
    tmp = canonical_form(i[1],sigma0)
    if tmp not in tree_dict:
        tree_dict[tmp] = i+[1]
    else:
        tree_dict[tmp][3]+= 1

cases = sorted(tuple(tree_dict[i] for i in tree_dict), key=lambda c: c[0])

print ("\nTotal cases",len(cases))
for i in range(len(cases)):
    print("Case",i,"maxdet = ",cases[i][0],"vertex of the additional simplex next to sigma0 = ",cases[i][1][-1].pts[0],"Chosen from",cases[i][3],"options")

print("\n")

#raise ValueError("STOP")


for i in cases:
    start = datetime.now()

    print("Working on vertex",i[1][-1].pts[0],".",i[3],"copies. maxdet =",i[0])


    # After we placed the simplex sigma next to sigma0, we need to make sure that the target is the proper type A neighbor
    # The expected free sum is known
    # 0, e4, (-1,0,0,0,1,1) 
    # e2, e3, (-1,-1,-1,1,1,1)
    # e5, e6, (-1,0,0,1,0,0) 
    # So we replace appropriate vertex
    if i[1][-1].pts[0]==(-1,0,0,0,1,1):
        target = i[1][-1].replace([0,0,0,0,0,1],[-1,0,0,1,0,0])
        target_index = i[1][-1].pts.index((0,0,0,0,0,1))
    elif i[1][-1].pts[0]==(-1,0,0,1,0,0):
        target = i[1][-1].replace([0,0,0,0,0,0],[-1,0,0,0,1,1])
        target_index = i[1][-1].pts.index((0,0,0,0,0,0))
    elif i[1][-1].pts[0]==(-1,-1,-1,1,1,1):
        target = i[1][-1].replace([0,0,0,0,0,1],[-1,0,0,1,0,0])
        target_index = i[1][-1].pts.index((0,0,0,0,0,1))
    else:
        raise ValueError("Unexpected vertex of sigma")

    target_facet = i[1][-1].facet(target_index)
    target_vertex = i[1][-1].pts[target_index]

    new_candidates = neighbor_candidates(i[1], target_facet, target_vertex, edges, bary)

    # removing target simplex from the candidates
    deviations = [i for i in new_candidates if symmetry_key(target.pts)!=symmetry_key(i.pts)]

    if len(new_candidates) - len(deviations) != 1:
        raise ValueError("The target simplex must be in the candidates")

    open_trees = []

    closed = {
    "frame":0,
    "open":0,
    "certified":0
    }


    # Checking all candidates and keeping only those that are not refuted by frame and not certified
    for j in deviations:
        tmp = process_child(i[1],j)
        closed[tmp[0]]+=1
        if tmp[0] == "open":
            tmp2=maxdet(tmp[1])[1]
            if tmp2 is None:
                raise ValueError("maxdet failed for a configuration")
            open_trees.append([tmp2,tmp[1],tmp[2]])

    if len(deviations) != sum(closed[j] for j in closed):
        raise ValueError("Every configuration must be processed")

    print("Total candidates other than target for this case:",len(deviations))
    print("Total closed by frame:",closed["frame"])
    print("Total certified:",closed["certified"])
    print("Total open remaining for this case:",closed["open"],"\n")

    # Creating dictionary of remaining configuration in this case up to symmetries keeping sigma0
    tree_dict = {}

    for j in open_trees:
        tmp = canonical_form(j[1],sigma0)
        if tmp not in tree_dict:
            tree_dict[tmp] = j+[1]
        else:
            tree_dict[tmp][3]+= 1

    new_cases = sorted(tuple(tree_dict[j] for j in tree_dict), key=lambda c: c[0])

    print ("\nTotal candidates in this case up to symmetry:",len(new_cases))
    print ("\n")

    for j in range(len(new_cases)):
    

        print("Working on candidate",j,"for this case.",new_cases[j][3],"copies. maxdet =",new_cases[j][0])


        # After we identified the target simplex of type A, we do not restrict other types of neighbors
        tmp = close_tree([[new_cases[j][1],new_cases[j][2]]],edges,bary,100)
    
        
        print("Total visited for candidate",j,"for this case:",tmp["total"])
        print("Total certified:",tmp["certified"])
        print("Total children:",tmp["children"])
        print("Total excluded by frame cutoff:",tmp["frame"])
        #print("Total excluded by library:",tmp["library"])
        print("Total duplicates:",tmp["duplicate"])
        #print("Total duplicate roots:",tmp["duplicate_roots"])
        print("Maximal configuration size:",tmp["max_size"])
        print("\n")
    end = datetime.now()
    print("Total time spent for this case: ", (end-start).total_seconds(),"\n")

print("All cases except target simplices have been discarded")
print("Conclusion: one free sum forces all neighboring free sums due to symmetries\n")

print("Total time:", (datetime.now()-total_start).total_seconds(),"seconds")


    
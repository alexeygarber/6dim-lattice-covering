# Setting
# We start from the standard simplex sigma0
# 
# Suppose all neighbors of sigma0 are of types B or C

# Claim ZZZ: then the determinant is below D0

# The only assumption is that all other neighbors are of types A, B, C


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



det = maxdet([sigma0])

if det[1] is None:
    raise ValueError("maxdet failed on the initial simplex")

if det[1]<float(D0):
    raise ValueError("This contradicts existence of Q_6^c1")


cases = [[det[1],[sigma0],det[3]]]




# For every facet of sigma0, attach a simplex of type B or C
for index in range(7):
    new_vertex = sigma0.pts[index]


    new_trees = []

    closed = {
    "frame":0,
    "open":0,
    "certified":0
    }

    children = 0

    for tree in cases:
        candidates = neighbor_candidates(tree[1], sigma0.facet(index), new_vertex, edges, bary_BC)

        children += len(candidates)

    # Checking all candidates and keeping only those that are not refuted by frame and not certified
        for i in candidates:
            tmp = process_child(tree[1],i)
            closed[tmp[0]]+=1
            if tmp[0] == "open":
                tmp2=maxdet(tmp[1])[1]
                if tmp2 is None:
                    raise ValueError("maxdet failed for a configuration")
                new_trees.append([tmp2,tmp[1],tmp[2],1,cases.index(tree)])

    if children != sum(closed[i] for i in closed):
        raise ValueError("Every configuration must be processed")

    print("Total candidates attached at facet",index,":",children)
    print("Total closed by frame:",closed["frame"])
    print("Total certified:",closed["certified"])
    print("Total open remaining:",closed["open"],"\n")

    # Creating dictionary of remaining configuration up to symmetry
    # We use new dictionary for each size of the tree
    tree_dict = {}

    for i in new_trees:
        # Here we are restricting symmetries to those that map sigma0 to sigma0
        tmp = canonical_form(i[1],sigma0)
        if tmp not in tree_dict:
            tree_dict[tmp] = i
        else:
            tree_dict[tmp][3]+= 1

    cases = sorted(tuple(tree_dict[i] for i in tree_dict), key=lambda c: c[0])

    print ("For face",index,"across",new_vertex,"Total cases",len(cases))
    for i in range(len(cases)):
        print("Case",i,"parent Case",cases[i][4],". maxdet = ",cases[i][0],"vertex = ",cases[i][1][-1].pts[0],"Chosen from",cases[i][3],"options")


    print("\n")

seen = set()

totals = {"total":0, 
            "children":0, 
            "frame":0, 
            "library":0, 
            "certified":0,
            "duplicate":0, 
            "duplicate_roots":0, 
            "max_size":0,
            "max_visited":0}

for i in cases:
    start = datetime.now()

    print("Working on case",cases.index(i),".",i[3],"copies. maxdet =",i[0])


    # After we placed all types B or C around sigma0, we do not restrict other types of neighbors
    tmp = close_tree([[i[1],i[2]]],edges,bary,100,seen=seen)
    
    end = datetime.now()
        
    print("Total visited for case",cases.index(i),":",tmp["total"])
    print("Total certified:",tmp["certified"])
    print("Total children:",tmp["children"])
    print("Total excluded by frame cutoff:",tmp["frame"])
    #print("Total excluded by library:",tmp["library"])
    print("Total duplicates:",tmp["duplicate"])
    #print("Total duplicate roots:",tmp["duplicate_roots"])
    print("Maximal configuration size:",tmp["max_size"])
    print("Total time spent: ", (end-start).total_seconds(),"\n")

    for key in totals:
        if key=="max_size":
            if tmp[key]>totals[key]:
                totals[key]=tmp[key]
        elif key=="max_visited":
            if tmp["total"]>totals[key]:
                totals[key]=tmp["total"]
        else:
            totals[key]+= tmp[key]


print("Total visited for all cases combined:",totals["total"])
print("Total certified:",totals["certified"])
print("Total children:",totals["children"])
print("Total excluded by frame cutoff:",totals["frame"])
#print("Total excluded by library:",totals["library"])
print("Total duplicates:",totals["duplicate"])
#print("Total duplicate roots:",totals["duplicate_roots"])
print("Maximal configuration size:",totals["max_size"])
print("Maximal visited in one case:",totals["max_visited"])

print("Total time:", (datetime.now()-total_start).total_seconds(),"seconds")


    
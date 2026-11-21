#
# Here we implement configuration of simplices
# In a valid configuration, each facet belongs to at most two simplices that are on opposite side of the facet
# Important: to avoid costly checks, we allow simplices to intersect "improperly" provided the facet condition holds and parity classes of edges do not clash
#

import simplex
import itertools
import numpy as np
from cutoff import EDGE_CUTOFF


# First an auxiliary function
# Computes barycentric coordinates x_i of the edge with respect to simplex and returns (sum x_i)^2 + sum(x_i^2)
def edge_coord_squared(S,edge):
    tmp = [sum([edge[i]*S.Vinv[i][j] for i in range(6)]) for j in range(6)]
    return sum(i*i for i in [sum(tmp)]+tmp)

# Checking if adding simplex T to a configuration of simplices violates EDGE_CUTOFF with respect oany of the new simplex placed as the standard one
def frame_violation(simplices, T):
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")
    if len({S.key for S in simplices}) != len(simplices):
        raise ValueError("Simplices must be different")
    if not isinstance(T,simplex.Simplex):
        raise ValueError("Only simplices allowed")
    for i in simplices:
        for j in T.edges:
            if edge_coord_squared(i,j)>EDGE_CUTOFF:
                return (i,j)
    for j in {e for i in simplices for e in i.edges}:
        if edge_coord_squared(T,j)>EDGE_CUTOFF:
            return (T,j)

    return None

# Checks that the configuration of simplices does not have two simplices on one side of each facet
# Note that we do not check whether simplices overlap
# Only opposite sides of one facet
def side_violation(simplices):
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")
    if len({S.key for S in simplices}) != len(simplices):
        raise ValueError("Simplices must be different")
    facets = {}
    for i in simplices:
        for j in range(7):
            tmp = simplex.symmetry_key(i.facet(j))
            if tmp in facets:
                if facets[tmp][0]>1:
                    raise ValueError("No facet may belong to more than two simplices of any configuration")

                facets[tmp][0] += 1
                # This means the common facets are translations
                if simplex.translation_key(i.facet(j))==simplex.translation_key(facets[tmp][1].facet(facets[tmp][2])):
                    tmp_pt = tuple(facets[tmp][1].pts[facets[tmp][2]][k]+min(i.facet(j))[k]-min(facets[tmp][1].facet(facets[tmp][2]))[k] for k in range(6))
                # Otherwise, they are reflections
                else:
                    tmp_pt = tuple(-facets[tmp][1].pts[facets[tmp][2]][k]+min(i.facet(j))[k]+max(facets[tmp][1].facet(facets[tmp][2]))[k] for k in range(6))
                    
                if i.bary(tmp_pt)[j]!=(-1):
                    return (facets[tmp][1],i)
            else:
                facets[tmp] = [1,i,j]

    return None


# From the list of types, create potential neighbors with respect to the standard simplex
# Types are A, B, C, D, E

def apex_coefficients(patterns):
    for i in patterns:
        if len(i)!=6:
            raise ValueError("Each pattern must have length 6")
        for j in i:
            if j!=int(j):
                raise ValueError("Each coordinate has to be an integer")
        if sum(i)!=2:
            raise ValueError("Sum of each pattern must be 2")

    tmp = []
    for i in patterns:
        tmp = tmp + [b for b in sorted(set(itertools.permutations(tuple(int(c) for c in i))))]
    
    return tmp

# dictionary of facets with their multiplicity in simplices
def facet_table(simplices):
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")
    if len({S.key for S in simplices}) != len(simplices):
        raise ValueError("Simplices must be different")
    facets = {}
    for i in simplices:
        for j in range(7):
            tmp = simplex.symmetry_key(i.facet(j))
            if tmp in facets:
                facets[tmp][0] += 1
            else:
                facets[tmp] = [1,i.facet(j),i.pts[j]]
    for i in facets:
        if facets[i][0]>2:
            raise ValueError("No facet may belong to more than two simplices of any configuration")
    return facets
    


# For the collection of simplices, gives all facets that belong to only one simplex
# All checks are done in facet_table
def open_facets(simplices):    
    facets = facet_table(simplices)
    return [tuple([facets[i][1],facets[i][2]]) for i in facets if facets[i][0]==1]




# Return a list of simplices that we will try to add to the configuration
# simplices = current simplices in the configuration
# facet = facet we are trying to close
# vertex = existing opposite vertex in the facet
# edges = potential edges
# bary = potential barycentric coordinates for a new vertex (the coordinate of vertex excluded)

def neighbor_candidates(simplices, facet, vertex, possible_edges, bary):
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")

    if len(facet)!=6:
        raise ValueError("Facet must have 6 vertices")
    
    for i in facet:
        if len(i)!=6:
            raise ValueError("Each vertex must have 6 coordinates")
        for j in i:
            if int(j)!=j:
                raise ValueError("Each coordinate must be an integer")

    if len(vertex)!=6:
        raise ValueError("Each vertex must have 6 coordinates")
    for j in vertex:
        if int(j)!=j:
            raise ValueError("Each coordinate must be an integer")    

    facet = sorted((tuple(int(j) for j in i) for i in facet))
    vertex = tuple(int(i) for i in vertex)

    if all([(set(facet)|{vertex})!=set(S.pts) for S in simplices]):
        raise ValueError("facet and vertex must give a simplex of the configuration")

    facets = facet_table(simplices)

    if facets[simplex.symmetry_key(facet)][0]!=1:
        raise ValueError("facet has to be open in simplices")

    for i in bary:
        if len(i)!=6:
            raise ValueError("We need barycentric coordinates for each vertex of the facet")
        for j in i:
            if int(j)!=j:
                raise ValueError("Barycentric coordinates must be integers")
        if sum(i)!=2:
            raise ValueError("Sum of barycentric coordinates must be 2")

    old_edges = tuple(edge for S in simplices for edge in S.edges)

    

    old_parity = {}
    for v in old_edges:
        if simplex.parity_class(v) not in old_parity:
            old_parity[simplex.parity_class(v)]=v


    keys = {S.key for S in simplices}

    out = []

    for ind in bary:
        new_vert = tuple(sum([ind[i]*facet[i][j] for i in range(6)])-vertex[j] for j in range(6))

        new_edges = tuple(simplex.sign_normalize(tuple(new_vert[j]-i[j] for j in range(6))) for i in facet)
        new_parity = tuple(simplex.parity_class(v) for v in new_edges)

        # Checking that all new edges are allowed by the cutoff
        
        if not all([i in possible_edges for i in new_edges]):
            continue

        # Checking that there is no parity collision between old and new edges
        if any([((new_parity[i] in old_parity) and new_edges[i]!= old_parity[new_parity[i]]) for i in range(6)]):
            continue



        # Checking if facet+vertex is a simplex in the existing configuration
        # If it is, then we treated closed facet as open
        new_simplex = simplex.Simplex([new_vert]+facet)
        if new_simplex.key in keys:
            raise ValueError("We cannot have another simplex containing a closed facet")

        # Checking that each new facet belongs to at most two simplices (one new and potentially one old)
        # i = 0 is not needed because facets(0) is the original facet
        tmp_key = [simplex.symmetry_key(new_simplex.facet(i)) for i in range(1,7)]

        if any((facets[i][0] if i in facets else 0)>1 for i in tmp_key):
            continue

        # Checking that if a facet of new simplex belongs to two simplices, then the simplices are on opposite sides
        tmp_sides = False
        for i in range(6):
            if (tmp_key[i] in facets) and (facets[tmp_key[i]][0]==1):
                # This check uses that simplices are unimodular
                # If two facets are translations of each other, we add vectors to apexes
                if simplex.translation_key(new_simplex.facet(i+1))==simplex.translation_key(facets[tmp_key[i]][1]):
                    tmp_vector = tuple( new_simplex.V[i][j] + facets[tmp_key[i]][2][j] - facets[tmp_key[i]][1][0][j] for j in range(6))
                # Otherwise, we subtract
                else:
                    tmp_vector = tuple( new_simplex.V[i][j] - facets[tmp_key[i]][2][j] + facets[tmp_key[i]][1][0][j] for j in range(6))
                # Then tmp_vector must lie in the subspace generated by the vectors of the common facet

                if simplex.b_det(
                    [tmp_vector]+list(new_simplex.V[:i])+list(new_simplex.V[(i+1):]))!=0:
                    tmp_sides=True

        if tmp_sides:
            continue


        out.append(new_simplex)


    return out

# Returns approximate circumradius squared of a potential new simplex with respect to the current approximate maximum Q
def circumradius_squared(new_simplex,Q):
    if not isinstance(new_simplex,simplex.Simplex):
        raise ValueError("Circumradius must be computed for a simples")
    Q = np.asarray(Q, dtype=float)
    if Q.shape!=(6,6):
        raise ValueError("Q must be 6x6")
    V = np.array(new_simplex.V, dtype=float)
    G = (V @ Q) @ V.transpose()
    diag = np.diag(G)
    return float(diag @ (np.linalg.solve(G, diag)))/4



# Choose an open facet from the collection of simplices
# Find all open facets
# For each facet, find potential simplices and their circumradius with respect to the currently best Q
# Choose the one with max min circumradius. If equal, with fewer candidates.

def choose_facet(simplices,Q,possible_edges,bary):

    if len(simplices)==0:
        raise ValueError("Cannot choose open facet for empty configuration")

    
    if Q is not None:
        Q = np.asarray(Q, dtype=float)
        if Q.shape!=(6,6):
            raise ValueError("Q must be 6x6 or None")

    candidates = {}

    for i in open_facets(simplices):
        candidates[i]=neighbor_candidates(simplices,i[0],i[1],possible_edges,bary)


    # length 0 means we don't have open facets
    if len(candidates)==0:
        # If we have exactly 360 simplices, then we found a complete triangulation.
        # Will need to optimize the covering radius over the L-type (if exists)
        if len(simplices)==360:
            raise ValueError("Complete triangulation found!")
        # Otherwise, we can discard
        else: 
            return (None,None,[])

    if Q is None:
        min_number = float("inf")
        for i in candidates:
            if len(candidates[i])<min_number:
                min_number=len(candidates[i])
                cand_facet = i
    else:
        min_R = 0
        min_number = float("inf")
        for i in candidates:
            tmp_R = min((circumradius_squared(j,Q) for j in candidates[i]),default=float("inf"))
            if (tmp_R > min_R+(1e-4)):
                min_R = tmp_R
                min_number = len(candidates[i])
                cand_facet = i
            # The second condition is needed to handle inf and inf
            if (abs(tmp_R - min_R)<=(1e-4)) or (tmp_R == min_R):
                if len(candidates[i])<min_number:
                    min_number=len(candidates[i])
                    cand_facet = i
    return (cand_facet[0],cand_facet[1],candidates[cand_facet])

if __name__ == "__main__":

    import edges
    import cutoff
    import relax

    L1=simplex.Simplex([[0,0,0,0,0,0],
                    [1,0,0,0,0,0],
                    [0,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])
    
    L2=simplex.Simplex([[2,1,1,0,-1,-1],
                    [1,0,0,0,0,0],
                    [0,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])

    L3=simplex.Simplex([[1,1,0,0,0,0],
                        [1,0,0,0,0,0],
                        [0,1,0,0,0,0],
                        [0,0,1,0,0,0],
                        [0,0,0,1,0,0],
                        [0,0,0,0,1,0],
                        [4,4,-1,-1,-1,1]])

    L4=simplex.Simplex([[1,1,0,0,0,0],
                            [1,0,0,0,0,0],
                            [0,1,0,0,0,0],
                            [0,0,1,0,0,0],
                            [0,0,0,1,0,0],
                            [0,0,0,0,1,0],
                            [0,0,0,0,0,1]])

    L5=simplex.Simplex([[1,1,1,1,0,-2],
                        [1,0,0,0,0,0],
                        [0,1,0,0,0,0],
                        [0,0,1,0,0,0],
                        [0,0,0,1,0,0],
                        [0,0,0,0,1,0],
                        [0,0,0,0,0,1]])
    

    root = [L1,L5]

    tmp = choose_facet(root,relax.maxdet(root)[3],edges.edge_candidates(cutoff.EDGE_CUTOFF),apex_coefficients([
 [1, 1, 0, 0, 0, 0],
 [1, 1, 1, 0, 0, -1],
 [1, 1, 1, 1, -1, -1],
 [1, 1, 1, 1, 0, -2],
 [2, 1, 1, 0, -1, -1]
 ]))

    print(len(tmp[2]))

    print(frame_violation([L1,L2,L3],L5))

    L6=simplex.Simplex([[1,1,-1,-1,0,0],
                    [-1,0,0,0,0,0],
                    [0,-1,0,0,0,0],
                    [0,0,-1,0,0,0],
                    [0,0,0,-1,0,0],
                    [0,0,0,0,-1,0],
                    [0,0,0,0,0,-1]])

    L7=simplex.Simplex([[0,0,0,0,0,2],
                    [1,0,0,0,0,2],
                    [0,1,0,0,0,2],
                    [1,1,-1,0,0,2],
                    [0,0,0,1,0,2],
                    [0,0,0,0,1,2],
                    [0,0,0,0,0,3]])

    # checking two cases when simplices are on one or on different sides
    print(side_violation([L1,L7]))
    print(side_violation([L1,L6]))
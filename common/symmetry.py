# Computes the canonical form of a collection of classes simplices
# Each class is up to translations and central reflections
# The canonical form does not depend on translations and GL_6 symmetries

import numpy as np
import itertools
import simplex

# All permutations on 7 elements
PERMS = np.array(list(itertools.permutations(range(7))))

# Powers of 128 for canonical form computation
POW = 128 ** np.arange(5, -1, -1)

#  We write difference between barycentric coordinates of the vertices with respect to the base simplex. So they sum to 0
#  These are coordinates that participate in our EDGE_CUTOFF
def bary_simplices(base,T):
    B = np.array([base.bary(p) for p in T.pts], dtype=np.int64)
    return B - B[0]

def canonical_form(simplices, base = None):
    if len(simplices)==0:
        raise ValueError("List of simplices must be nonempty")
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Input needs to be a list of simplices")
    if base is not None:
        if not isinstance(base,simplex.Simplex):
            raise ValueError("Base has to be a simplex")
        if not(base.key in {S.key for S in simplices}):
            raise ValueError("Base needs to be in the list of simplices")
    if len(simplices)!= len({S.key for S in simplices}):
        raise ValueError("All simplices in the list must be different")
    
    if base is None:
        bases = simplices
    else:
        bases = [base]


    # For each base S, this is list of barycentric coordinates of edges with respect to S
    tmp = [[bary_simplices(S,T) for T in simplices] for S in bases]

    global_min = None

    for i in range(len(bases)):
        # Choose only the coordinates corresponding to bases[i]        
        D = np.array(tmp[i])

        # Verify that our coordinates are not too large. This is used to make sure that sorting of integers below coincides with the lexicographic order. The bound is guaranteed by the EDGE_CUTOFF bound
        if np.abs(D).max() >= 32:
            raise ValueError("Edges cannot have length more than 31. This violates the EDGE_CUTOFF bound")

        # For given bases Construct an array indexed by
        # simplex in simplices (7x7 barycentric coordinates of edges) and permutations (7! potential reordering of vertices of the base simplex)
        # img[a,b,c,d] is constructed as follows
        # c = ordering of the base simplex
        # d = position in that ordering (with the initial removed)
        # a = simplex
        # b = vertex of simplex a ( minus the initial vertex of simplex a)
        img = D[:, :, PERMS[:, 1:]]

        # We use powers of 128 to "sort" coordinate representations. Since each entry is less than 32 in absolute value, their products with the vector POW sort according to the lexicographic sort of vectors. For each of 7! choices of the order of vertices from bases[i], we get 7 numbers corresponding to barycentric coordinates of vertex - vertex0 of each simplex T
        # Properties we need
        # 1) Bijection (because entries of D are small)
        # 2) Order is the same as lexicographic order
        # 3) code (x-y) = code(x) - code(y): preservation of subtraction
        codes = img @ POW
        
        # For translative key
        # For each ordering of the verties of the base simplex, we choose minimal vertex of the simplex and subtract it from all (similarly to our translative key in simplex.py)
        # Then we sort the vertices for each simplex and each ordering. In the end we get translative key of a simplex for each ordering
        # 
        # Similarly for the reflection key
        # 
        tr_key = np.sort(codes - codes.min(axis=1, keepdims=True), axis=1)
        ref_key = np.sort(codes.max(axis=1, keepdims=True) - codes, axis=1)


       
        diff = tr_key - ref_key
        # Choose the indices where the keys are different
        first = np.argmax(diff != 0, axis=1)

        # Pick first nonzero entries in diff
        d = np.take_along_axis(diff, first[:, None, :], axis=1)[:, 0, :]

        # Pick those that are smaller
        # As the outcome, for each order of the vertices of the base, key contains symmetry_key of each simplex in simplices and these keys are in the same order as simplices
        key = np.where((d < 0)[:, None, :], tr_key, ref_key)
        
        # Changing the order of indices to
        # 1) Ordering
        # 2) simplex
        # 3) key
        key = key.transpose(2,0,1)

        # For each ordering of vertices, sorting the keys lexicographically and keeping the minimal key only
        min_key = min(tuple(sorted(map(tuple, k))) for k in key.tolist())


        # comparing the minimum key for this basis simplex to the current best minimum
        if global_min is None:
            global_min = min_key
        else:
            if min_key < global_min:
                global_min = min_key


        


    return global_min

if __name__ == "__main__":

    #import edges
    #import cutoff
    #import relax

    L1=simplex.Simplex([[0,0,0,0,0,0],
                    [1,0,0,0,0,0],
                    [0,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])

    L2=simplex.Simplex([[1,1,0,0,0,0],
                    [1,0,0,0,0,0],
                    [0,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])    
    
    L3=simplex.Simplex([[2,2,2,1,-2,-3],
                        [1,0,0,0,0,0],
                        [0,1,0,0,0,0],
                        [0,0,1,0,0,0],
                        [0,0,0,1,0,0],
                        [0,0,0,0,1,0],
                        [0,0,0,0,0,1]])    

    print(canonical_form([L1,L2,L3]))

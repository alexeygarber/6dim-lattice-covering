import simplex
import symmetry
from configuration import frame_violation, choose_facet, side_violation
from certificate import certificate_bound
from relax import maxdet, rational_psd
from cutoff import D0
from datetime import datetime


# We assume simplices are frame-checked
def certify_child(simplices):
    if len(simplices)==0:
        raise ValueError("The list of simplices should be nonempty")
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")
    if len(simplices)!= len({S.key for S in simplices}):
        raise ValueError("All simplices in the list must be different")

    child_maxdet = maxdet(simplices)
    if (child_maxdet[1] is not None) and (child_maxdet[1] < float(D0)):
        rational_matrices = [rational_psd(M,10**6) for M in child_maxdet[2]]
        child_cert = certificate_bound(simplices, rational_matrices)
        if (child_cert is not None) and (child_cert < D0):
            support = [(S, R) for S, R in zip(simplices, rational_matrices) if any(i!=0 for j in R for i in j)]
            
            # Checking that removing non-supporting simplices does not change the certified bound
            if len(support)<len(simplices):
                if certificate_bound([i[0] for i in support],[i[1] for i in support])!=child_cert:
                    raise ValueError("Removing zero matrices cannot change the certified bound")
            return (True, child_maxdet[3], support)

    return (False, child_maxdet[3], None)


# This is the order of the checks
# First check if we get too long edges
# Then we check if a subset is in the library
# And third is maxdet bound + certification


def process_child(parent_simplices,new_simplex):
    # checks happen in frame_violation
    check_frame = frame_violation(parent_simplices,new_simplex)
    if check_frame is not None:
        return ("frame",check_frame)

    child = list(parent_simplices)+[new_simplex]

    # Will check for library here
    #
    #
    #

    cert = certify_child(child)
    if cert[0]:
        return ("certified",cert[2])

    return ("open", child, cert[1])
    




# Default for report_every is 100
# If a noninteger or negative is used, we report nothing
def close_tree(root_configurations,edges,bary, report_every = 100, duplicate_check=True, seen = None):
    if not (isinstance(report_every, int) and report_every > 0):
        report_every = 0

    if seen is None:
        seen = set()

    all_stack = list(root_configurations)
    
    for i in all_stack:
        if frame_violation_all(i[0]) is not None:
            raise ValueError("Each root configuration must not violate edge lengths")
        if side_violation(i[0]) is not None:
            raise ValueError("Each root configuration must not have two simplices sharing a facet from one side")
        if simplex.parity_collision(tuple(edge for S in i[0] for edge in S.edges)) is not None:            
            raise ValueError("There is a parity collision in an initial configuration")

    stack = []
        
    for i in all_stack:
        tmp = symmetry.canonical_form(i[0])
        if tmp not in seen:
            stack.append(i)
            seen.add(tmp)

    different = len(stack)

    duplicate_roots = len(all_stack) - different
    
    total_visited = 0
    total_frame_closed = 0
    total_library = 0
    total_certified = 0
    total_duplicate = 0
    total_children = 0

    max_size = 0
    start = datetime.now()    
    while len(stack)>0:
        

        if report_every != 0:
            if total_visited% report_every==0:
                current = datetime.now()
                print("Total visited:",total_visited, ". Stack size: ",len(stack),". Working on these roots: ",(current-start).total_seconds(),"seconds")

        node = stack.pop()
        total_visited +=1

        if len(node[0])> max_size:
            max_size = len(node[0])

        chosen_facet = choose_facet(node[0], node[1], edges, bary)
        for i in chosen_facet[2]:
            total_children += 1

            child_outcome = process_child(node[0],i)

            if child_outcome[0]=="frame":                
                total_frame_closed +=1

            elif child_outcome[0]=="certified":                
                total_certified +=1

            elif child_outcome[0]=="library":                
                total_library +=1

            elif child_outcome[0] == "open":
                if duplicate_check:
                    tmp = symmetry.canonical_form(child_outcome[1])
                    if tmp in seen:
                        total_duplicate += 1
                    else:
                        seen.add(tmp)
                        stack.append([child_outcome[1],child_outcome[2]])
                else:
                    stack.append([child_outcome[1],child_outcome[2]])

            else:
                raise ValueError("child_outcome produced unexpected result")


    if total_children - total_certified - total_frame_closed - total_library - total_duplicate != total_visited - different:
        raise ValueError("The tree data does not match")



    return {"total":total_visited, 
            "children":total_children, 
            "frame":total_frame_closed, 
            "library":total_library, 
            "certified":total_certified,
            "duplicate":total_duplicate, 
            "duplicate_roots":duplicate_roots, 
            "max_size":max_size}








def frame_violation_all(simplices):
    if len(simplices)==0:
        raise ValueError("The list of simplices should be nonempty")
    for i in simplices:
        if not isinstance(i,simplex.Simplex):
            raise ValueError("Only simplices allowed")
    if len(simplices)!= len({S.key for S in simplices}):
        raise ValueError("All simplices in the list must be different")
    for k in range(1,len(simplices)):
        tmp = frame_violation(simplices[:k],simplices[k])
        if tmp is not None:
            return tmp
    return None

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

    print(certify_child([L1,L3]))
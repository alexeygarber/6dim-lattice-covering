import itertools
import math

from simplex import parity_collision
from cutoff import EDGE_CUTOFF

def edge_candidates(R):
    """
    For given R, creates the list of all edges (one of plus-minus pair is chosen) that may theoretically appear in a Delaunay triangulation with the standard simplex and with sum of squares of entries and square of sum of entries at most R
    """
    if R<=0:
        raise ValueError("Expected bound should be positive")
    tmp = math.isqrt(R)

    edge_set=[
        (1,0,0,0,0,0),(0,1,0,0,0,0),(0,0,1,0,0,0),(0,0,0,1,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1),
        (1,-1,0,0,0,0),(1,0,-1,0,0,0),(1,0,0,-1,0,0),(1,0,0,0,-1,0),(1,0,0,0,0,-1),
        (0,1,-1,0,0,0),(0,1,0,-1,0,0),(0,1,0,0,-1,0),(0,1,0,0,0,-1),
        (0,0,1,-1,0,0),(0,0,1,0,-1,0),(0,0,1,0,0,-1),
        (0,0,0,1,-1,0),(0,0,0,1,0,-1),
        (0,0,0,0,1,-1),
    ]

    return frozenset(tuple(p) for p in itertools.product(range(-tmp,tmp+1), repeat=6) if ((sum([i*i for i in p])+(sum([i for i in p]))**2)<=R and tuple(p)>(tuple(-i for i in p)) and parity_collision(edge_set+[tuple(p)]) is None))

# Testing


if __name__ == "__main__":
    
    print(len(edge_candidates(EDGE_CUTOFF)))
from functools import cached_property

def b_det(matrix):
    """
    Computes the exact integer determinant of a square matrix 
    using the Bareiss fraction-free algorithm.
    """
    # Create a copy 
    M = [list(row) for row in matrix]
    n = len(M)
    
    det_sign = 1
    prev_pivot = 1
    
    for k in range(n):
        # 1. Row Pivoting: If the diagonal element is 0, find a non-zero row below it
        if M[k][k] == 0:
            pivot_row = -1
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    pivot_row = r
                    break
            
            # If no non-zero pivot is found, the determinant is 0
            if pivot_row == -1:
                return 0
                
            # Swap rows and flip the sign of the determinant
            M[k], M[pivot_row] = M[pivot_row], M[k]
            det_sign *= -1
        
        pivot = M[k][k]
        
        # 2. Update the remaining submatrix elements
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                # Bareiss update formula: (M[i][j]*M[k][k] - M[i][k]*M[k][j]) / prev_pivot
                # Using integer division (//) since remainder is guaranteed to be 0

                tmp = divmod(M[i][j] * pivot - M[i][k] * M[k][j], prev_pivot)

                if tmp[1]!=0:
                    raise ValueError("Bareiss det failed. Remainder is not 0")

                M[i][j] = tmp[0]
        
        # Keep track of the current pivot for the next iteration's denominator
        prev_pivot = pivot

    # The bottom-right entry of the modified upper-triangular matrix is the determinant
    return int(det_sign * M[n - 1][n - 1])
        
def sign_normalize(v):
    if all([int(i)==0 for i in v]):
        raise ValueError("Zero vector cannot be sign normalized")

    return max(tuple((int(i) for i in v)),tuple((int(-i) for i in v)))

def parity_class(v):
    if len(v)!=6:
        raise ValueError("Parity classes are considered for vectors with 6 entries only")
    return tuple(int(i%2) for i in v)

def parity_collision(vecs):
    """
    The function checks if there are parity collisions in a list/tuple of vectors.

    If there is an even vector with parity class (0,0,0,0,0,0), then it returns this vector (provided it is not zero)
    
    If there are two vectors u and v in the same parity class that are not equal or opposite, then it returns the first pair.

    It raises error through sign_normalize if there is a zero vector in the list

    Otherwise, it returns None
    """

    vecs = tuple(tuple(int(j) for j in i) for i in vecs)

    if not (all(len(i)==6 for i in vecs)):
        raise ValueError("A vector submitted for parity check does not have length 6")

    fixed = {}

    for v in vecs:
        tmp_parity = parity_class(v)
        if tmp_parity == (0,0,0,0,0,0):
            return (sign_normalize(v),)
        if tmp_parity not in fixed:
            fixed[tmp_parity]=sign_normalize(v)
        else:
            if fixed[tmp_parity]!=sign_normalize(v):
                return (fixed[tmp_parity],sign_normalize(v))
    return None


def translation_key(points):
    points = tuple(tuple(c for c in p) for p in points)
    for p in points:
        if len(p)!=6:
            raise ValueError("Each point must have six coordinates")
        for c in p:
            if int(c) != c:
                raise ValueError("Each point must have integer coordinates")

    tmp = tuple(tuple(int(c) for c in p) for p in points)

    min_pt = min(tmp)

    return tuple(sorted(tuple(i[j]-min_pt[j] for j in range(6)) for i in tmp))

def symmetry_key(points):
    points = tuple(tuple(c for c in p) for p in points)
    for p in points:
        if len(p)!=6:
            raise ValueError("Each point must have six coordinates")
        for c in p:
            if int(c) != c:
                raise ValueError("Each point must have integer coordinates")
    
    return min(translation_key(points),translation_key(tuple(tuple(-int(c) for c in p) for p in points)))

    

class Simplex:
    def __init__(self, pts):
        for p in pts:
            for c in p:
                if int(c) != c:
                    raise ValueError("Each simplex must have integer vertices")

        self.pts = tuple(tuple(int(c) for c in p) for p in pts)
        if len(pts)!=7:
            raise ValueError("Simplex must have 7 vertices")
        for p in pts:
            if len(p)!=6:
                raise ValueError("Each vertex of simplex must have 6 entries")
        
        self.V = tuple(tuple(self.pts[i][j]-self.pts[0][j] for j in range(6)) for i in range(1,7))
        self.detV = b_det(self.V)

        if abs(self.detV)!=1:
            raise ValueError("Only unimodular simplices are allowed")


    # Computing inverse matrix
    @cached_property
    def Vinv(self):
        return tuple(
            tuple(self.detV*((-1)**(i+j))*b_det([self.V[k][:j]+self.V[k][j+1:] for k in range(6) if k!=i]) for i in range(6)) for j in range(6)
        )

    def facet(self, j):
        return (self.pts[:j] + self.pts[j+1:])

    def bary(self,x):
    
        if len(x)!=6:
            raise ValueError("Point has to have six coordinates")
        x=tuple((int(i) for i in x))
        tmp = tuple(sum([(x[i]-self.pts[0][i])*self.Vinv[i][j] for i in range(6)]) for j in range(6))
        
        return (1-sum(tmp),)+tuple(tmp)

    def point(self,w):

        if len(w)!=7:
            raise ValueError("There must be 7 baricentric coordinates")

        w= tuple((int(i) for i in w))

        if sum(w)!=1:
            raise ValueError("Barycentric coordinates should add to 1")

        return tuple([sum([self.pts[j][i]*w[j] for j in range(7)]) for i in range(6)])


    # Returning a neighbor simplex by replacing vertex v with vertex w
    def replace(self, v, w):
        if len(v)!=6:
            raise ValueError("Each point must have 6 coordinates")
        if len(w)!=6:
            raise ValueError("Each point must have 6 coordinates")
        for c in v:
            if int(c) != c:
                raise ValueError("Each point must have integer coordinates")
        for c in w:
            if int(c) != c:
                raise ValueError("Each point must have integer coordinates")  
        v = tuple(int(c) for c in v)
        if v not in self.pts:
            raise ValueError("Replace point must be a vertex of the simplex")

        i = self.pts.index(v)
        w = tuple(int(c) for c in w)

        if self.bary(w)[i]!=(-1):
            raise ValueError("w should be on the other side of v w.r.t. common facet")

        return Simplex(list(self.pts[:i])+[w]+list(self.pts[i+1:]))                      
        

    @cached_property
    def edges(self):
        return tuple(sign_normalize(tuple(self.pts[j][k]-self.pts[i][k] for k in range(6))) for i in range(7) for j in range(i+1,7))

    @cached_property
    def key(self):
        return symmetry_key(self.pts)    

    #__eq__, __hash__

if __name__ == "__main__":

    #import edges
    #import cutoff
    #import relax

    L1=Simplex([[0,0,0,0,0,0],
                    [1,0,0,0,0,0],
                    [0,1,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,1,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1]])

    
    

    print(L1.replace([0,0,0,0,0,0],[1,1,0,0,0,0]).pts)

    print(L1.replace([0,1,0,0,0,0],[1,1,0,0,0,0]).pts)
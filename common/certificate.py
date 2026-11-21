# We are checking that we given restrictions, determinants of all quadratic forms are at most D0=4.395397
# Goal: pi^3/(6*det(Lambda)) \leq 2.4649
# This implies: det(Q) \geq pi^6/(6*2.4649)^2 > 4.3953971
# 

from fractions import Fraction as F
import simplex
import math

def certificate_bound(simplices,Ms):
    if len(simplices)==0:
        raise ValueError("List of simplices has to be nonempty")
        
    for i in simplices:
        if not isinstance(i, simplex.Simplex):
            raise ValueError("Only Simplex allowed in certificate bound")

    if len(Ms)!=len(simplices):
        raise ValueError("Each simplex must have assigned matrix")

    for i in Ms:
        if len(i) != 7:
            raise ValueError("Each assigned matrix has to be 7-by-7")
        for j in i:
            if len(j) != 7:
                raise ValueError("Each assigned matrix has to be 7-by-7")
            for k in j:
                if not isinstance(k,F):
                    raise ValueError("Entries of each matrix M have to be fractions")

        # Checking if this M is PSD
        
        tmp = [[0]*7 for _ in range(7)]

        for k in range(7):
            for j in range(7):
                if i[j][k]!=i[k][j]:
                    raise ValueError("Each matrix M has to be symmetric")
                tmp[j][k]=i[j][k]

        for k in range(7):
            if tmp[k][k]<0:
                raise ValueError("Each matrix M has to be PSD")
            if tmp[k][k]==0:
                for j in range(k+1,7):
                    if tmp[j][k]!=0:
                        raise ValueError("Each matrix M has to be PSD")
            if tmp[k][k]>0:
                for j in range(k+1,7):
                    for l in range(k+1,7):
                        tmp[j][l]-=tmp[j][k]*tmp[k][l]/tmp[k][k]
        
    

    C = sum([4*i[0][0] for i in Ms])

    # This is the right-hand side
    if C<=0:
        return None

    tmp = [[0]*6 for _ in range(6)]

    for i in range(len(Ms)):
        tmp2 = [[0]*6 for _ in range(6)]

        # Creating matrix M'+2D
        for k in range(6):
            tmp2[k][k] += 2*Ms[i][0][k+1]
            for l in range(6):
                tmp2[k][l] += Ms[i][k+1][l+1]
    
        
        # Multiplying V (M'+2D) V^t and adding (actually subtracting) together for all simplices. Note that we store vertices as rows.
        for k in range(6):
            for l in range(6):
                tmp[k][l] -= sum([sum([simplices[i].V[a][k]*tmp2[a][b]*simplices[i].V[b][l] for a in range(6)]) for b in range(6)])        

    den = math.lcm(*(x.denominator for row in tmp for x in row))

    for i in range(6):
        for j in range(6):
            tmp[i][j]=int(den*tmp[i][j])

    tmp = tuple([tuple(i) for i in tmp])

    for j in range(5):
        if simplex.b_det(tuple(i[:(j+1)] for i in tmp[:(j+1)]))<=0:
            return None

    det = simplex.b_det(tmp)

    if det<=0:
        return None
    else:    
        return (C/6)**6/F(det,(den)**6)

if __name__ == "__main__":

    import relax
    
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

    tmp = relax.maxdet([L1,L2])
    psd = [relax.rational_psd(i,10**2) for i in tmp[2]]

    print(tmp[1])

    print(tmp[2])

    print(psd)

    tmp_bound = certificate_bound([L1,L2],psd)

    print(float(tmp_bound))
    print(tmp_bound)


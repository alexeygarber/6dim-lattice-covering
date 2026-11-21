from fractions import Fraction as F
#from mpmath import mp


# This is our target density
Theta0=F(24649,10000)

# We are searching for PQF with determinant at least D0
D0=F(4395397,1000000)


# Verification using approximation of pi
PI_LOWER = F(3141592653589793238462643383279, 10**30)

if D0>(((PI_LOWER)**3)/(6*Theta0))**2:
    raise ValueError("Cutoffs loaded incorrectly")

# Cutoff of potential edges
# (x_1+...+x_6)^2+x_1^2+...+x_6^2 <= 31
EDGE_CUTOFF = 31
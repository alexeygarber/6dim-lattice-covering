# This script computes maxdet for all potential neighbors of the standard simplex
# The remaining simplices are of five types
# B = (1,1,0,0,0,0)
# A = (1,1,1,0,0,-1)
# C = (1,1,1,1,-1,-1)
# E = (1,1,1,1,0,-2)
# D = (2,1,1,0,-1,-1)


from pathlib import Path
import sys
import os
import json



sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))

from simplex import Simplex
import relax
import certificate
from cutoff import D0

script_dir = Path(__file__).resolve().parent
with open(script_dir / "step2_apexes.json") as f:
    
    apexes = [[int(x) for x in a] for a in json.load(f)]

for i in apexes:
    if len(i)!=6:
        raise ValueError("Each potential apex should have 6 coordinates")
    if sum(i)!=2:
        raise ValueError("Sum of coordinates of each potential apex must be 2")

sigma0=Simplex([
    [0,0,0,0,0,0],
    [1,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,1,0,0],
    [0,0,0,0,1,0],
    [0,0,0,0,0,1]
])

remaining = []  

file_path = script_dir / "step2_ABCDE.json"
file_path.unlink(missing_ok=True)

for i in apexes:
    sigma1=Simplex([
    i,
    [1,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,1,0,0],
    [0,0,0,0,1,0],
    [0,0,0,0,0,1]
    ])

    tmp = relax.maxdet([sigma0,sigma1])

    outcome = True
    cert = None

    if tmp[1] is not None:
        if tmp[1]<float(D0):
            cert = certificate.certificate_bound([sigma0,sigma1],[relax.rational_psd(M,10**6) for M in tmp[2]])
            if (cert is not None) and (cert<D0):
                outcome = False

    if outcome:
        remaining.append(i)
        print(i,"- Accepted. Solver status =", tmp[0]+".","Convex optimization bound=",tmp[1])
    else:
        print(i,"- Rejected. Solver status =", tmp[0]+".","Certified bound=",float(cert))
    


with open(file_path, "w") as f:
    f.write("[\n" + ",\n".join(json.dumps(t) for t in remaining) + "\n]\n")

with open(file_path, "r") as f:
    data = json.load(f)
    if ([t for t in data]!=remaining):
        raise ValueError("Writing to file was not successful")

print("Total remaining",len(remaining),"candidates")
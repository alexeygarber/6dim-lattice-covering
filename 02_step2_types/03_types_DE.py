from pathlib import Path
import sys
import os
import json
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common'))

from simplex import Simplex
from relax import maxdet
from cutoff import EDGE_CUTOFF
from edges import edge_candidates
from configuration import apex_coefficients
from search import close_tree

script_dir = Path(__file__).resolve().parent
with open(script_dir / "step2_ABCDE.json") as f:
    
    apexes = [[int(x) for x in a] for a in json.load(f)]

for i in apexes:
    if len(i)!=6:
        raise ValueError("Each potential apex should have 6 coordinates")
    if sum(i)!=2:
        raise ValueError("Sum of coordinates of each potential apex must be 2")

if [2,1,1,0,-1,-1] not in apexes:
    raise ValueError("Type D must be in the apexes")

if [1,1,1,1,0,-2] not in apexes:
    raise ValueError("Type E must be in the apexes")

edges = edge_candidates(EDGE_CUTOFF)
bary = apex_coefficients(apexes)


file_path_remaining = script_dir / "step2_ABC.json"
file_path_remaining.unlink(missing_ok=True)
    
for i in [[2,1,1,0,-1,-1],[1,1,1,1,0,-2]]:
    start = datetime.now()
    #tmp = search_DE(i)

    print("Working on type",i)

    sigma0=Simplex([
    [0,0,0,0,0,0],
    [1,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,1,0,0],
    [0,0,0,0,1,0],
    [0,0,0,0,0,1]])

    sigma1=Simplex([
    i,
    [1,0,0,0,0,0],
    [0,1,0,0,0,0],
    [0,0,1,0,0,0],
    [0,0,0,1,0,0],
    [0,0,0,0,1,0],
    [0,0,0,0,0,1]])

    tmp = close_tree([[[sigma0,sigma1],maxdet([sigma0,sigma1])[3]]],edges,bary, report_every = 10, duplicate_check = False, seen = None)


    end = datetime.now()
    print("Total visited for",i,":",tmp["total"],"\nTotal certified:",tmp["certified"],"\nTotal children:",tmp["children"],"\nTotal excluded by frame cutoff:",tmp["frame"],"\nMaximal configuration size:",tmp["max_size"],"\nTotal time spent: ", (end-start).total_seconds(),"\n")
    apexes.remove(i)
        

with open(file_path_remaining, "w") as f:
    f.write("[\n" + ",\n".join(json.dumps(t) for t in apexes) + "\n]\n")

with open(file_path_remaining, "r") as f:
    data = json.load(f)
    if ([t for t in data]!=apexes):
        raise ValueError("Writing to file was not successful")
import json, glob, numpy as np
from collections import defaultdict
size=defaultdict(list); powr=defaultdict(list)
for f in sorted(glob.glob("calib_result_*.json")):
    d=json.load(open(f))
    for k,v in d["size"].items(): size[k].append(v)
    for k,v in d["power_fine"].items(): powr[k].append(v)
print("=== SIZE (debiased ref), mean over seeds, target 0.05 ===")
for k in sorted(size,key=float):
    a=np.array(size[k]); print(f"  ws={k}: {a.mean():.4f} (se {a.std()/len(a)**.5:.4f})")
print("=== POWER FINE grid at ws=0.75, mean over seeds ===")
for k in sorted(powr,key=lambda x:float(x[2:])):
    a=np.array(powr[k]); print(f"  {k}: {a.mean():.4f} (se {a.std()/len(a)**.5:.4f})")

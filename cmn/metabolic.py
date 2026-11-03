"""Task 2.3 part 3 -- ATP-weighted metabolic network.
Usage: python metabolic.py [path/to/metabolic_energy_balance_network.csv]"""
import csv, sys
from shortest_paths import Graph, dijkstra, bellman_ford, path_to, cycle_weight, INF

PATH = sys.argv[1] if len(sys.argv) > 1 else "data/metabolic_energy_balance_network.csv"
SOURCE = "Glucose"
TARGETS = ["Pyruvate", "Oxaloacetate", "Serine"]

def load(path):
    names, ids, edges = [], {}, []
    def vid(x):
        if x not in ids:
            ids[x] = len(names); names.append(x)
        return ids[x]
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            edges.append((vid(r["Source"].strip()), vid(r["Target"].strip()), int(r["ATP Cost"])))
    g = Graph(len(names), directed=True)            # reactions are directed
    for u, v, w in edges:
        g.add_edge(u, v, w)
    return g, names, ids

def fmt(d):
    return "-inf" if d == -INF else ("unreachable" if d == INF else f"{d:+d}")

def simple_path_costs(g, s):
    """Brute force over SIMPLE paths (fine for 25 edges): min cost per vertex."""
    best = {s: 0}
    def dfs(u, cost, seen):
        for v, w in zip(g.adj[u], g.wt[u]):
            if v in seen: continue
            if v not in best or cost + w < best[v]: best[v] = cost + w
            dfs(v, cost + w, seen | {v})
    dfs(s, 0, {s})
    return best

g, names, ids = load(PATH)
s = ids[SOURCE]
print(f"{g.n} metabolites, {sum(len(a) for a in g.adj)} reactions; weights in "
      f"{sorted({w for ws in g.wt for w in ws})}\n")

# ---- 1. Dijkstra (assumes non-negative weights)
dd, pd = dijkstra(g, s)
print("DIJKSTRA from", SOURCE)
for t in TARGETS:
    print(f"  {t:13s} distance {fmt(dd[ids[t]]):>4s}   path: " + " -> ".join(names[x] for x in path_to(pd, s, ids[t])))

# ---- 2. Bellman-Ford
db, pb, cyc = bellman_ford(g, s)
print("\nBELLMAN-FORD from", SOURCE)
for t in TARGETS:
    print(f"  {t:13s} distance {fmt(db[ids[t]])}")
print("\n  negative-weight cycle reachable from the source:", cyc is not None)
if cyc:
    print("  cycle:", " -> ".join(names[x] for x in cyc))
    for a, b in zip(cyc, cyc[1:]):
        w = [w for v, w in zip(g.adj[a], g.wt[a]) if v == b][0]
        print(f"      {names[a]:>24s} -> {names[b]:<24s} {w:+d}")
    print(f"  total weight of the cycle: {cycle_weight(g, cyc):+d} ATP per lap")

# ---- 3. why Dijkstra is wrong, concretely
print("\nWHY DIJKSTRA CANNOT BE TRUSTED HERE")
pep, pyr = ids["Phosphoenolpyruvate"], ids["Pyruvate"]
print(f"  Dijkstra settles Phosphoenolpyruvate at {dd[pep]:+d} (path "
      + " -> ".join(names[x] for x in path_to(pd, s, pep)) + ").")
lap = path_to(pd, s, ids["Oxaloacetate"])
print("  But it then reaches Oxaloacetate at", fmt(dd[ids["Oxaloacetate"]]), "and the edge Oxaloacetate -> PEP costs +1,")
print(f"  so PEP is reachable at {dd[ids['Oxaloacetate']] + 1:+d} < {dd[pep]:+d}: a settled vertex gets a cheaper path, which Dijkstra's invariant forbids.")
print("  Going round the cycle again lowers every vertex downstream of PEP by 1 ATP per lap, without bound.")

# ---- 4. full table: Bellman-Ford vs Dijkstra vs best simple path
sp = simple_path_costs(g, s)
print("\nALL METABOLITES  (distance from Glucose)")
print(f"  {'metabolite':28s} {'Dijkstra':>9s} {'Bellman-Ford':>13s} {'best simple path':>17s}")
for v in range(g.n):
    flag = "  <-- Dijkstra wrong" if dd[v] != db[v] else ""
    print(f"  {names[v]:28s} {fmt(dd[v]):>9s} {fmt(db[v]):>13s} {fmt(sp.get(v, INF)):>17s}{flag}")
wrong = [names[v] for v in range(g.n) if dd[v] != db[v]]
print(f"\nDijkstra disagrees with Bellman-Ford on {len(wrong)} of {g.n} metabolites.")

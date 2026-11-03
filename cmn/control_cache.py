"""Control experiment: same n and m as regime A, but neighbours are i+-1, i+-2
(sequential memory access) instead of random.  If cache misses explain the
super-linear slope of BFS on random graphs, time per vertex stays flat here."""
import csv, gc, random, time
from shortest_paths import Graph, bfs, dijkstra

def lattice(n):
    g = Graph(n)
    for i in range(n):
        g.add_edge(i, (i + 1) % n, 1); g.add_edge(i, (i + 2) % n, 1)
    return g

def randomg(n, rng):
    g = Graph(n)
    for v in range(1, n): g.add_edge(rng.randrange(v), v, 1)
    for _ in range(n + 1):
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v: g.add_edge(u, v, 1)
    return g

def tm(f, g):
    best = 1e9
    for _ in range(2):
        gc.collect(); gc.disable(); t = time.perf_counter(); f(g, 0); best = min(best, time.perf_counter() - t); gc.enable()
    return best

rows = []
print(f"{'n':>8} | {'BFS random us/vertex':>21} | {'BFS lattice us/vertex':>22}")
for n in (10**3, 10**4, 10**5, 10**6):
    tr = tm(bfs, randomg(n, random.Random(1))); tl = tm(bfs, lattice(n))
    rows.append(dict(n=n, us_per_vertex_random=tr / n * 1e6, us_per_vertex_lattice=tl / n * 1e6))
    print(f"{n:8d} | {tr/n*1e6:21.3f} | {tl/n*1e6:22.3f}", flush=True)
with open("results/control_cache.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

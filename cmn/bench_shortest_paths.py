"""Scaling benchmark for BFS / Dijkstra.  Usage: python bench_shortest_paths.py
Regime A: mean degree 4, n = 1e3 .. 1e6.   Regime B: density 0.01, n = 1e2 .. 1e4.
Writes results/shortest_paths.csv"""
import csv, gc, os, random, time
from heapq import heappush, heappop
from shortest_paths import Graph, bfs, dijkstra, INF

SEED = 424242
WMAX = 100

def graph_mean_degree(n, deg, rng):
    """Random recursive tree (guarantees connectivity) + random extra edges until
    m = deg*n/2 undirected edges.  Parallel edges are possible, self-loops are not."""
    g = Graph(n)
    for v in range(1, n):
        g.add_edge(rng.randrange(v), v, rng.randint(1, WMAX))
    for _ in range(deg * n // 2 - (n - 1)):
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v:
            g.add_edge(u, v, rng.randint(1, WMAX))
    return g

def graph_density(n, dens, rng):
    """G(n, m) with m = dens * n(n-1)/2 distinct edges (no forced connectivity)."""
    m = round(dens * n * (n - 1) / 2)
    g, seen = Graph(n), set()
    while len(seen) < m:
        u, v = rng.randrange(n), rng.randrange(n)
        if u == v: continue
        key = (u, v) if u < v else (v, u)
        if key not in seen:
            seen.add(key); g.add_edge(u, v, rng.randint(1, WMAX))
    return g

def dijkstra_counting(g, s):
    """Same as dijkstra() but counts heap pushes (not timed)."""
    dist = [INF] * g.n; done = [False] * g.n; dist[s] = 0
    heap = [(0, s)]; pushes = 1
    while heap:
        d, u = heappop(heap)
        if done[u]: continue
        done[u] = True
        for v, w in zip(g.adj[u], g.wt[u]):
            nd = d + w
            if nd < dist[v] and not done[v]:
                dist[v] = nd; heappush(heap, (nd, v)); pushes += 1
    return pushes

def best_time(f, g, s, reps):
    best = 1e18
    for _ in range(reps):
        gc.collect(); gc.disable()
        t = time.perf_counter(); f(g, s); best = min(best, time.perf_counter() - t)
        gc.enable()
    return best

def run(regime, n, make, rows):
    rng = random.Random(SEED + n)
    t0 = time.perf_counter(); g = make(n, rng); tb = time.perf_counter() - t0
    m_arcs = sum(len(a) for a in g.adj)
    # source: highest-degree vertex (lies in the giant component when there is one)
    s = max(range(n), key=lambda u: len(g.adj[u]))
    reps = 3 if n <= 10**5 else 2
    t_bfs = best_time(bfs, g, s, reps)
    t_dij = best_time(dijkstra, g, s, reps)
    dist, _ = bfs(g, s)
    reached = sum(d < INF for d in dist)
    pushes = dijkstra_counting(g, s)
    rows.append(dict(regime=regime, n=n, m=m_arcs // 2, reached=reached, t_bfs=t_bfs,
                     t_dijkstra=t_dij, heap_pushes=pushes))
    print(rows[-1], f"(build {tb:.1f}s)", flush=True)

def main():
    os.makedirs("results", exist_ok=True)
    rows = []
    for n in (10**3, 3 * 10**3, 10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6):
        run("A: mean degree 4", n, lambda n, r: graph_mean_degree(n, 4, r), rows)
    for n in (100, 200, 500, 1000, 2000, 5000, 10000):
        run("B: density 0.01", n, lambda n, r: graph_density(n, 0.01, r), rows)
    with open("results/shortest_paths.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

if __name__ == "__main__":
    main()

"""Validate BFS / Dijkstra / Bellman-Ford against networkx on random graphs.
networkx is used ONLY here, to check answers.  Usage: python validate_shortest_paths.py"""
import random
import networkx as nx
from shortest_paths import Graph, bfs, dijkstra, bellman_ford, path_to, cycle_weight, INF

SEED = 20261102

def random_graph(rng, directed, wlo, whi):
    n = rng.randint(2, 150)
    p = rng.choice([0.01, 0.02, 0.05, 0.1, 0.3, 0.6])      # varying density
    g, G = Graph(n, directed), (nx.DiGraph() if directed else nx.Graph())
    G.add_nodes_from(range(n))
    for u in range(n):
        # directed: all ordered pairs; undirected: each unordered pair once
        for v in (range(n) if directed else range(u + 1, n)):
            if u != v and rng.random() < p:
                w = rng.randint(wlo, whi)
                g.add_edge(u, v, w); G.add_edge(u, v, weight=w)
    return g, G, rng.randrange(n)

def check_paths(g, s, dist, parent, weighted):
    """Every returned path must be a real path whose cost equals dist[v]."""
    bad = 0
    for v in range(g.n):
        if dist[v] == INF:
            bad += path_to(parent, s, v) is not None
            continue
        p = path_to(parent, s, v)
        if p is None or p[0] != s or p[-1] != v:
            bad += 1; continue
        cost = 0
        for a, b in zip(p, p[1:]):
            ws = [w for x, w in zip(g.adj[a], g.wt[a]) if x == b]
            if not ws: bad += 1; break
            cost += min(ws) if weighted else 1
        else:
            bad += cost != dist[v]
    return bad

def main():
    rng = random.Random(SEED)
    dis_bfs = dis_dij = dis_bf = bad_paths = 0
    for t in range(100):
        g, G, s = random_graph(rng, directed=(t % 2 == 1), wlo=1, whi=20)
        ref_h = nx.single_source_shortest_path_length(G, s)
        ref_d = nx.single_source_dijkstra_path_length(G, s)
        dh, ph = bfs(g, s); dd, pd = dijkstra(g, s)
        dis_bfs += any(dh[v] != ref_h.get(v, INF) for v in range(g.n))
        dis_dij += any(dd[v] != ref_d.get(v, INF) for v in range(g.n))
        bad_paths += check_paths(g, s, dh, ph, False) + check_paths(g, s, dd, pd, True)
    print(f"BFS      : {dis_bfs} disagreements on 100 graphs")
    print(f"Dijkstra : {dis_dij} disagreements on 100 graphs")
    print(f"invalid paths returned (BFS+Dijkstra): {bad_paths}")

    # Bellman-Ford: negative weights, with and without reachable negative cycles
    cyc_found = cyc_nx = 0
    for t in range(100):
        g, G, s = random_graph(rng, directed=True, wlo=-5, whi=15)
        db, pb, cyc = bellman_ford(g, s)
        try:
            ref = nx.single_source_bellman_ford_path_length(G, s)
            nxcyc = False
        except nx.NetworkXUnbounded:
            nxcyc = True
        cyc_nx += nxcyc; cyc_found += cyc is not None
        if (cyc is not None) != nxcyc:
            dis_bf += 1; continue
        if cyc is None:
            dis_bf += any(db[v] != ref.get(v, INF) for v in range(g.n))
        else:
            dis_bf += cyc[0] != cyc[-1] or cycle_weight(g, cyc) >= 0
            # vertices NOT marked -inf must have the exact distance: recompute with
            # networkx on the subgraph without the marked vertices (it has no
            # reachable negative cycle, and no shortest path to an unmarked
            # vertex can pass through a marked one)
            keep = [v for v in range(g.n) if db[v] != -INF]
            if s in keep:
                ref2 = nx.single_source_bellman_ford_path_length(G.subgraph(keep), s)
                dis_bf += any(db[v] != ref2.get(v, INF) for v in keep)
            else:   # s is on/after the cycle: everything reachable from s must be -inf
                dis_bf += any(db[v] not in (-INF, INF) for v in range(g.n))
    print(f"Bellman-Ford: {dis_bf} disagreements on 100 graphs "
          f"({cyc_nx} had a reachable negative cycle, mine flagged {cyc_found}; "
          f"every printed cycle verified to be closed and negative)")

if __name__ == "__main__":
    main()

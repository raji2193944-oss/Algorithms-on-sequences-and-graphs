import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shortest_paths import Graph, bfs, dijkstra, bellman_ford, path_to, cycle_weight, INF

def small():
    g = Graph(6, directed=True)
    for u, v, w in [(0,1,7),(0,2,9),(0,5,14),(1,2,10),(1,3,15),(2,3,11),(2,5,2),(3,4,6),(5,4,9)]:
        g.add_edge(u, v, w)
    return g

def test_dijkstra_textbook():
    d, p = dijkstra(small(), 0)
    assert d == [0, 7, 9, 20, 20, 11]
    assert path_to(p, 0, 4) == [0, 2, 5, 4]

def test_bfs_hops_ignore_weights():
    d, p = bfs(small(), 0)
    assert d == [0, 1, 1, 2, 2, 1]

def test_unreachable_and_source():
    g = Graph(4, directed=True); g.add_edge(0, 1, 3)
    for f in (bfs, dijkstra):
        d, p = f(g, 0)
        assert d[2] == INF and d[3] == INF and path_to(p, 0, 3) is None
        assert path_to(p, 0, 0) == [0]
    assert bellman_ford(g, 0)[0][1] == 3

def test_zero_weight_edges_and_ties():
    g = Graph(4, directed=True)
    for u, v, w in [(0,1,0),(1,2,0),(0,2,0),(2,3,5),(1,3,5)]:
        g.add_edge(u, v, w)
    assert dijkstra(g, 0)[0] == [0, 0, 0, 5]

def test_undirected_edge_stored_both_ways():
    g = Graph(3); g.add_edge(0, 1, 4); g.add_edge(1, 2, 1)
    assert dijkstra(g, 2)[0] == [5, 1, 0]

def test_bellman_ford_negative_edges_no_cycle():
    g = Graph(4, directed=True)
    for u, v, w in [(0,1,4),(0,2,5),(1,2,-3),(2,3,2)]:
        g.add_edge(u, v, w)
    d, p, c = bellman_ford(g, 0)
    assert d == [0, 4, 1, 3] and c is None

def test_bellman_ford_negative_cycle():
    g = Graph(5, directed=True)
    for u, v, w in [(0,1,1),(1,2,-1),(2,3,-1),(3,1,1),(0,4,2)]:   # cycle 1->2->3->1 = -1
        g.add_edge(u, v, w)
    d, p, c = bellman_ford(g, 0)
    assert c is not None and c[0] == c[-1] and set(c) == {1, 2, 3}
    assert cycle_weight(g, c) == -1
    assert d[0] == 0 and d[4] == 2 and d[1] == d[2] == d[3] == -INF

def test_dijkstra_terminates_on_negative_cycle():
    g = Graph(3, directed=True)
    for u, v, w in [(0,1,1),(1,2,-1),(2,1,-1)]:
        g.add_edge(u, v, w)
    dijkstra(g, 0)          # must not loop forever (wrong answer is expected, hanging is not)

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok  ", name)

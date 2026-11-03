"""Task 2.3 -- BFS, Dijkstra, Bellman-Ford on adjacency lists.

Representation: vertices are 0..n-1.  adj[u] is the list of neighbours of u and
wt[u] the parallel list of edge weights (wt[u][i] is the weight of u->adj[u][i]).
An undirected edge is stored twice.  Memory O(n + m).
"""
import math
from collections import deque
from heapq import heappush, heappop

INF = math.inf


class Graph:
    def __init__(self, n, directed=False):
        self.n = n
        self.directed = directed
        self.adj = [[] for _ in range(n)]
        self.wt = [[] for _ in range(n)]

    def add_edge(self, u, v, w=1):
        self.adj[u].append(v)
        self.wt[u].append(w)
        if not self.directed:
            self.adj[v].append(u)
            self.wt[v].append(w)

    def edges(self):
        for u in range(self.n):
            for v, w in zip(self.adj[u], self.wt[u]):
                yield u, v, w


def path_to(parent, s, t):
    """Vertex list s..t following parent pointers, or None if t was not reached."""
    if t != s and parent[t] == -1:
        return None
    p = [t]
    while p[-1] != s:
        p.append(parent[p[-1]])
    return p[::-1]


def bfs(g, s):
    """Hop distances from s (weights ignored).  O(n + m): each vertex is
    enqueued at most once and each adjacency list is scanned once."""
    dist = [INF] * g.n
    parent = [-1] * g.n
    dist[s] = 0
    q = deque([s])
    adj = g.adj
    while q:
        u = q.popleft()
        du = dist[u] + 1
        for v in adj[u]:
            if dist[v] == INF:
                dist[v] = du
                parent[v] = u
                q.append(v)
    return dist, parent


def dijkstra(g, s):
    """Textbook Dijkstra with a binary heap (heapq) and lazy deletion.
    Correct for NON-NEGATIVE weights.  A vertex is settled the first time it is
    popped and is never touched again; with negative weights it still
    terminates, but its answer is not meaningful (see metabolic.py)."""
    n = g.n
    dist = [INF] * n
    parent = [-1] * n
    done = [False] * n
    dist[s] = 0
    heap = [(0, s)]
    adj, wt = g.adj, g.wt
    while heap:
        d, u = heappop(heap)
        if done[u]:                       # stale entry (lazy deletion)
            continue
        done[u] = True
        for v, w in zip(adj[u], wt[u]):
            nd = d + w
            if nd < dist[v] and not done[v]:
                dist[v] = nd
                parent[v] = u
                heappush(heap, (nd, v))
    return dist, parent


def bellman_ford(g, s):
    """Bellman-Ford, O(n*m).  Returns (dist, parent, cycle).

    dist[v] is the exact shortest-WALK distance from s; it is -inf for every
    vertex reachable from a negative cycle that is itself reachable from s
    (for those, no shortest walk exists).  cycle is None if no negative cycle
    is reachable from s, otherwise one such cycle as a vertex list
    [c0, c1, ..., c0] with edges c_i -> c_{i+1}."""
    n = g.n
    adj, wt = g.adj, g.wt
    dist = [INF] * n
    parent = [-1] * n
    dist[s] = 0
    for _ in range(n - 1):
        changed = False
        for u in range(n):
            du = dist[u]
            if du == INF:
                continue
            for v, w in zip(adj[u], wt[u]):
                if du + w < dist[v]:
                    dist[v] = du + w
                    parent[v] = u
                    changed = True
        if not changed:                   # converged early: no reachable negative cycle
            return dist, parent, None

    # still relaxable after n-1 passes  <=>  a negative cycle is reachable
    bad = [v for u in range(n) if dist[u] < INF
           for v, w in zip(adj[u], wt[u]) if dist[u] + w < dist[v]]
    if not bad:
        return dist, parent, None

    # (1) every vertex reachable from a relaxable vertex has distance -inf
    seen = set(bad)
    stack = list(bad)
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    final = list(dist)
    for v in seen:
        final[v] = -INF

    # (2) extract one cycle: do pass n on copies, walk parent pointers n steps
    d2, p2 = list(dist), list(parent)
    x = -1
    for u in range(n):
        if d2[u] == INF:
            continue
        for v, w in zip(adj[u], wt[u]):
            if d2[u] + w < d2[v]:
                d2[v] = d2[u] + w
                p2[v] = u
                x = v
    for _ in range(n):
        x = p2[x]
    cycle = [x]
    y = p2[x]
    while y != x:
        cycle.append(y)
        y = p2[y]
    cycle.append(x)
    cycle.reverse()
    return final, parent, cycle


def cycle_weight(g, cycle):
    """Total weight of the closed walk (uses the cheapest parallel edge)."""
    total = 0
    for a, b in zip(cycle, cycle[1:]):
        ws = [w for v, w in zip(g.adj[a], g.wt[a]) if v == b]
        if not ws:
            raise ValueError(f"no edge {a}->{b}")
        total += min(ws)
    return total

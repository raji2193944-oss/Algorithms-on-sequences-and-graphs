"""Task 2.3 part 3: Dijkstra vs Bellman-Ford on the ATP metabolic network.

Loads data/metabolic_energy_balance_network.csv (Source, Target, ATP Cost),
builds a directed Graph over integer vertices with a name<->index map, runs
Dijkstra from Glucose (numbers not meaningful here), then Bellman-Ford, and
reports the distances to Pyruvate, Oxaloacetate and Serine, plus a reachable
negative cycle if one exists.  Writes results/metabolic.txt.

All file I/O is explicit UTF-8: the CSV contains non-ASCII metabolite names
(α-Ketoglutarate), which Windows' default cp1252 mangles on read and cannot
encode at all on write.
"""
import csv
from shortest_paths import (Graph, dijkstra, bellman_ford, INF)

CSV_PATH = "data/metabolic_energy_balance_network.csv"
SOURCE = "Glucose"
TARGETS = ["Pyruvate", "Oxaloacetate", "Serine"]


def load(path=CSV_PATH):
    edges = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            edges.append((r["Source"], r["Target"], int(r["ATP Cost"])))
    names = sorted({u for u, _, _ in edges} | {v for _, v, _ in edges})
    idx = {name: i for i, name in enumerate(names)}
    g = Graph(len(names), directed=True)
    for u, v, w in edges:
        g.add_edge(idx[u], idx[v], w)
    return g, names, idx


def fmt(d):
    if d == -INF:
        return "-inf (reachable from a negative cycle)"
    if d == INF:
        return "unreachable"
    return str(int(d))


def cycle_weight(g, cycle):
    total = 0
    for a, b in zip(cycle, cycle[1:]):
        ws = [w for v, w in zip(g.adj[a], g.wt[a]) if v == b]
        if not ws:
            raise ValueError(f"no edge {a}->{b}")
        total += min(ws)
    return total


def main():
    g, names, idx = load()
    s = idx[SOURCE]
    out = []

    dd, _ = dijkstra(g, s)
    out.append("=== Dijkstra from Glucose (NOT trustworthy here) ===")
    for t in TARGETS:
        out.append(f"  {t}: {fmt(dd[idx[t]])}")
    out.append("  (three negative edges exist, so Dijkstra's assumption fails)")

    bf, _, cyc = bellman_ford(g, s)
    out.append("")
    out.append("=== Bellman-Ford from Glucose (correct) ===")
    for t in TARGETS:
        out.append(f"  {t}: {fmt(bf[idx[t]])}")
    if cyc is None:
        out.append("  no negative cycle reachable from Glucose")
    else:
        out.append(f"  reachable negative cycle (weight {cycle_weight(g, cyc)}):")
        out.append("    " + " -> ".join(names[v] for v in cyc))

    text = "\n".join(out)
    print(text)
    with open("results/metabolic.txt", "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()

# Algorithms-on-sequences-and-graphs
# CMN Coding Assignment 1 – algorithms on sequences and graphs

Python 3.9+ (developed on 3.13). The algorithmic core of every task uses the
standard library only (`heapq`, `collections.deque`). `networkx` is used **only**
in `validate_shortest_paths.py` to check answers. All random seeds are fixed in
the scripts (listed in each file's header).

The metabolic CSV contains non-ASCII metabolite names (α-Ketoglutarate); the
scripts open it with explicit `encoding="utf-8"`.

**Run every command from the `cmn/` directory of this repository** (the folder
containing the scripts, `data/` and `results/`). All paths below are relative
to `cmn/`.

```
pip install -r requirements.txt      # matplotlib, numpy, networkx (pytest optional)
```

## How to run the tests

```
python run_tests.py                  # everything, ~30 s including the networkx validation
```

or one file at a time (each also works with `pytest tests/`):

```
python tests/test_common_substring.py
python tests/test_palindromes.py
python tests/test_shortest_paths.py
python validate_shortest_paths.py    # 100 undirected/directed graphs vs networkx, plus
                                     # 100 Bellman-Ford graphs; prints the number of
                                     # disagreements (expected: 0 everywhere)
```

## Layout

| Task | Implementation | Tests | Benchmark → plot |
|---|---|---|---|
| 2.1 common k-mer | `common_substring.py` (`naive`, `hashing`) | `tests/test_common_substring.py` | `bench_common_substring.py`, `crossover_common_substring.py` → `plot_common_substring.py` |
| 2.2 revcomp palindromes | `palindromes.py` (`naive_pairs`, `expand_pairs`, `expand`) | `tests/test_palindromes.py` | `bench_palindromes.py` → `plot_palindromes.py` |
| 2.3 shortest paths | `shortest_paths.py` (`bfs`, `dijkstra`, `bellman_ford`) | `tests/test_shortest_paths.py`, `validate_shortest_paths.py` | `bench_shortest_paths.py`, `control_cache.py` → `plot_shortest_paths.py` |
| 2.3 part 3 | `metabolic.py` (reads `data/metabolic_energy_balance_network.csv`) | (no benchmark; single fixed network) | |

Outputs (CSV + PNG) go to `results/`.

## How to run the benchmarks and reproduce the results

```
# Task 2.1  (~4 min; naive is time-capped at 3 s per run and extrapolated as n^2 beyond that,
#            dashed in the plot -- pure-Python naive at n = 10^6 would take about a day, ~25-27 h)
python bench_common_substring.py
python crossover_common_substring.py       # fine scan of the crossover at small n
python plot_common_substring.py

# Task 2.2  -- needs the E. coli FASTA; if absent, fetch it:
#   curl -L -o data/u00096.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=fasta&retmode=text"
python bench_palindromes.py data/u00096.fa
python plot_palindromes.py
#   Refuses to run without a FASTA; use --synthetic for the control genome only.
#   Each run analyses the real genome AND a synthetic twin (same length and GC
#   content, generated from the i.i.d. model of the expected-count argument)
#   and prints the obs/exp comparison at the end.  The results shipped in this
#   archive are from the real E. coli chromosome (U00096.3); the prefix data
#   are prefixes of the real genome.

# Task 2.3  (~1.5 min; the largest graph, n = 10^6, needs ~1 GB RAM)
python validate_shortest_paths.py
python bench_shortest_paths.py
python control_cache.py
python plot_shortest_paths.py
python metabolic.py                        # Dijkstra vs Bellman-Ford on the ATP network
```

Timings are wall-clock `time.perf_counter()` on one machine; expect different absolute
numbers on yours (the slopes and crossovers are what matter).

## Conventions worth knowing before reading the code

* **2.1** both functions return `(i, j)` with `x[i:i+k] == y[j:j+k]`, or `None`.
* **2.2** a reverse-complement palindrome always has even length; the output "pairs" are all
  `(i, l)` with `l >= k` and `t[i:i+l]` a palindrome; "positions" is the set of starts `i`.
  Any character outside `ACGT` (e.g. `N`) never matches. Input is upper-cased.
* **2.3** graphs are `Graph(n, directed)` with parallel lists `adj[u]` / `wt[u]`;
  an undirected edge is stored twice. `dijkstra` settles each vertex once and is only
  correct for non-negative weights (on negative weights it terminates with a wrong
  answer; `bellman_ford` is the one to use, and returns `-inf` for vertices reachable
  from a negative cycle plus one such cycle).

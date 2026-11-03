""""Benchmark for Task 2.1.  Usage: python bench_common_substring.py

Writes results/common_substring.csv (the plot script makes the PNGs).
Rows are written to disk the moment they are measured, so an interrupted
run still leaves a partial (plottable) CSV.
"""
import csv, random, statistics, time, os
from common_substring import naive, hashing

SEED = 12345
AA = "ACDEFGHIKLMNPQRSTVWY"
SIZES = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10**4, 10**5, 10**6]
KS = [3, 10, 50]
CAP = 3.0            # stop measuring naive in a series once one run exceeds this (s)
REPS = 3

def rand_seq(rng, n, alphabet=AA):
    return "".join(rng.choices(alphabet, k=n))

def instance(family, n, k, rng):
    if family == "random":
        return rand_seq(rng, n), rand_seq(rng, n)
    # adversarial for naive: x = A^n, y = (A^{k-1} C)^*, no common k-mer,
    # every window comparison runs O(k) steps before it fails
    x = "A" * n
    y = (("A" * (k - 1) + "C") * (n // k + 1))[:n]
    return x, y

def timeit(f, x, y, k):
    t = time.perf_counter(); f(x, y, k); return time.perf_counter() - t

def main():
    os.makedirs("results", exist_ok=True)
    fieldnames = ["family", "k", "n", "t_hash", "t_naive"]
    with open("results/common_substring.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for family in ("random", "adversarial"):
            for k in KS:
                naive_alive = True
                for n in SIZES:
                    if n < k:
                        continue
                    tn, th = [], []
                    for r in range(REPS):
                        rng = random.Random(f"{SEED}-{r}-{n}-{k}")   # tuple: no collisions
                        x, y = instance(family, n, k, rng)
                        th.append(timeit(hashing, x, y, k))
                        if naive_alive:
                            tn.append(timeit(naive, x, y, k))
                        if family == "adversarial":
                            break    # deterministic instance, 1 rep is enough
                    row = dict(family=family, k=k, n=n,
                               t_hash=statistics.median(th),
                               t_naive=statistics.median(tn) if tn else "")
                    if tn and row["t_naive"] > CAP:
                        naive_alive = False
                    w.writerow(row)
                    f.flush()
                    print(row, flush=True)

if __name__ == "__main__":
    main()
### run this to produce the CSV, then run cmn/plot_common_substring.py to produce the PNGs
##or run "C:\Users\Sandra\AppData\Local\Programs\Python\Python313\python.exe cmn\plot_common_substring.py" in the terminal to produce the PNGs directly from the CSV
"""Fine scan for the crossover on random inputs. Uses many repetitions because
runs at n < 200 take microseconds. Writes results/crossover.csv"""
import csv, random, time
from common_substring import naive, hashing
AA = "ACDEFGHIKLMNPQRSTVWY"

def best_time(f, x, y, k, loops):
    best = 1e9
    for _ in range(5):
        t = time.perf_counter()
        for _ in range(loops): f(x, y, k)
        best = min(best, (time.perf_counter() - t) / loops)
    return best

rows = []
for k in (10, 50):
    for n in sorted(set(list(range(k, k + 6)) + [k + 10, k + 20, 2*k, 3*k, 5*k, 100, 150, 200])):
        rng = random.Random(7 + n + k)
        # average over several instances so one lucky instance doesn't decide
        tn = th = 0
        inst = 5
        for _ in range(inst):
            x = "".join(rng.choices(AA, k=n)); y = "".join(rng.choices(AA, k=n))
            tn += best_time(naive, x, y, k, 50); th += best_time(hashing, x, y, k, 50)
        rows.append(dict(k=k, n=n, t_naive=tn/inst, t_hash=th/inst, winner="naive" if tn < th else "hash"))
        print(rows[-1])
with open("results/crossover.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

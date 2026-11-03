import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open("results/common_substring.csv")))
fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
for a, fam in enumerate(("random", "adversarial")):
    for b, k in enumerate((3, 10, 50)):
        ax = axs[a][b]
        R = [r for r in rows if r["family"] == fam and int(r["k"]) == k]
        ns = [int(r["n"]) for r in R]
        ax.loglog(ns, [float(r["t_hash"]) for r in R], "o-", label="hashing")
        Rn = [r for r in R if r["t_naive"] != ""]
        ax.loglog([int(r["n"]) for r in Rn], [float(r["t_naive"]) for r in Rn], "s-", label="naive")
        # dashed quadratic extrapolation where naive was cut off by the time cap
        if len(Rn) < len(R) and not (fam == "random" and k == 3):
            n0, t0 = int(Rn[-1]["n"]), float(Rn[-1]["t_naive"])
            ext = [n for n in ns if n >= n0]
            ax.loglog(ext, [t0 * (n / n0) ** 2 for n in ext], "s--", color="C1", alpha=.5,
                      label="naive (extrapolated ~n²)")
        ax.set_title(f"{fam}, k={k}"); ax.grid(alpha=.3, which="both")
        if a == 1: ax.set_xlabel("n = m")
        if b == 0: ax.set_ylabel("time (s)")
        ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig("results/common_substring.png", dpi=150)

# empirical exponents from the measured points
import math
print("empirical log-log slopes between measured points (n>=1000)")
for fam in ("random", "adversarial"):
    for k in (3, 10, 50):
        R = [r for r in rows if r["family"] == fam and int(r["k"]) == k and int(r["n"]) >= 1000]
        for alg in ("t_hash", "t_naive"):
            pts = [(int(r["n"]), float(r[alg])) for r in R if r[alg] != ""]
            if len(pts) >= 2:
                s = math.log(pts[-1][1] / pts[0][1]) / math.log(pts[-1][0] / pts[0][0])
                print(f"  {fam:11s} k={k:2d} {alg:8s} slope {s:.2f}  (n {pts[0][0]}..{pts[-1][0]})")

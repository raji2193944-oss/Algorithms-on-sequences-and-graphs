import csv, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = list(csv.DictReader(open("results/palindromes_genome.csv")))
P = list(csv.DictReader(open("results/palindromes_prefix.csv")))
W = list(csv.DictReader(open("results/palindromes_worst.csv")))

slope = lambda x, y: math.log(y[-1] / y[0]) / math.log(x[-1] / x[0])
N_GENOME = 4_641_652

real = [r for r in G if r["genome"] == "real"]
syn  = [r for r in G if r["genome"] == "synthetic"]
ks = [int(r["k"]) for r in real]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))

# panel 0: counts vs k, real vs synthetic vs iid expectation
ax[0].semilogy(ks, [float(r["centres"]) for r in real], "o-",
               label="real genome (E. coli K-12)")
ax[0].semilogy(ks, [float(r["centres"]) for r in syn], "s-",
               label="synthetic twin (iid, same GC)")
ax[0].semilogy(ks, [float(r["expected_centres"]) for r in real], "^--",
               color="grey", label="expected, (n-k+1)·q^(k/2)")
ax[0].set_xlabel("k"); ax[0].set_ylabel("qualifying centres")
ax[0].set_xticks(ks)
ax[0].set_title("Palindromes vs k: real genome vs random baseline")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")

# panel 1: prefix scaling
for k in (4, 6, 8, 10):
    R = [r for r in P if int(r["k"]) == k]
    ns = [int(r["n"]) for r in R]
    ax[1].loglog(ns, [float(r["t_naive"]) for r in R], "s-", label=f"naive k={k}", alpha=.8)
R = [r for r in P if int(r["k"]) == 4]
ax[1].loglog([int(r["n"]) for r in R], [float(r["t_expand"]) for r in R],
             "o-", color="k", label="expansion k=4")
ax[1].set_xlabel("prefix length n"); ax[1].set_ylabel("time (s)")
ax[1].set_title("Real-genome prefixes: naive ~ n², expansion ~ n")
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")

# panel 2: worst-case families
for fam, c in (("A^(n/2)T^(n/2)", "C0"), ("(AT)^(n/2)", "C1")):
    R = [r for r in W if r["family"] == fam]
    ns = [int(r["n"]) for r in R]
    ax[2].loglog(ns, [float(r["t_naive"]) for r in R], "s-", color=c, label=f"naive, {fam}")
    ax[2].loglog(ns, [float(r["t_expand"]) for r in R], "o--", color=c, label=f"expansion, {fam}")
ax[2].set_xlabel("n"); ax[2].set_title("Worst-case families (k=4)")
ax[2].legend(fontsize=7); ax[2].grid(alpha=.3, which="both")

plt.tight_layout(); plt.savefig("results/palindromes.png", dpi=150)

# ---- console summary
print("prefix agreement: enforced by assertion in the bench (run would have "
      f"crashed otherwise); {len(P)} prefix runs all matched")
for k in (4, 6, 8, 10):
    R = [r for r in P if int(r["k"]) == k]
    ns = [int(r["n"]) for r in R]
    sn = slope(ns, [float(r["t_naive"]) for r in R])
    se = slope(ns, [float(r["t_expand"]) for r in R])
    print(f"k={k:2d} slope naive {sn:.2f}  expansion {se:.2f}")

R = [r for r in P if int(r["k"]) == 4 and int(r["n"]) == 8000][0]
est = float(R["t_naive"]) * (N_GENOME / 8000) ** 2
print(f"naive extrapolated to full genome (k=4): {est:.3g} s = {est/86400:.1f} days")

for fam in ("A^(n/2)T^(n/2)", "(AT)^(n/2)"):
    R = [r for r in W if r["family"] == fam]
    ns = [int(r["n"]) for r in R]
    print(fam, "slope naive %.2f  expansion %.2f  (pairs %s)" % (
        slope(ns, [float(r["t_naive"]) for r in R]),
        slope(ns, [float(r["t_expand"]) for r in R]),
        [int(r["pairs"]) for r in R]))

print("real genome, time per k:", [round(float(r["seconds"]), 2) for r in real])
print("real/synthetic obs_exp:",
      [f"{float(r['obs_over_exp']):.3f}" for r in real],
      "vs", [f"{float(r['obs_over_exp']):.3f}" for r in syn])

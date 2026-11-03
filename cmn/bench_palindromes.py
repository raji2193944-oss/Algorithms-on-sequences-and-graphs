"""Benchmark for Task 2.2 -- reverse-complement palindromes.

Usage:
    python bench_palindromes.py data/u00096.fa    # real genome + synthetic twin
    python bench_palindromes.py --synthetic       # synthetic genome only

Writes results/palindromes_genome.csv, _prefix.csv and _worst.csv.
"""
import csv, os, random, sys, time

from palindromes import naive_pairs, expand, summarize_pairs

SEED = 2026
KS = [4, 6, 8, 10]
PREFIX_SIZES = (500, 1000, 2000, 4000, 8000)
WORST_SIZES = (50, 100, 200, 400, 800)


def load_fasta(path):
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def gc_of(g):
    c = {b: g.count(b) for b in "ACGT"}
    return (c["C"] + c["G"]) / sum(c.values())


def synthetic(n, gc):
    rng = random.Random(SEED)
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]   # A C G T, strand-symmetric
    return "".join(rng.choices("ACGT", weights=w, k=n))


def t_of(f, *a):
    t0 = time.perf_counter()
    r = f(*a)
    return time.perf_counter() - t0, r


def save(name, rows):
    with open(f"results/{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    os.makedirs("results", exist_ok=True)

    # ---- pick the genomes to analyse
    if len(sys.argv) > 1 and not sys.argv[1].startswith("--"):
        g = load_fasta(sys.argv[1])
        genomes = [("real", g, sys.argv[1]),
                   ("synthetic", synthetic(len(g), gc_of(g)),
                    "synthetic twin (same length and GC)")]
    elif "--synthetic" in sys.argv:
        genomes = [("synthetic", synthetic(4_641_652, 0.508),
                    "synthetic random genome")]
    else:
        sys.exit("no FASTA given.\n"
                 "  python bench_palindromes.py data/u00096.fa\n"
                 "  python bench_palindromes.py --synthetic")

    # ---- (a) expansion on every genome; expected counts from its own composition
    rows = []
    for label, g, src in genomes:
        n = len(g)
        c = {b: g.count(b) for b in "ACGT"}
        tot = sum(c.values())
        pA, pC, pG, pT = (c[b] / tot for b in "ACGT")
        q = pA * pT * 2 + pC * pG * 2          # P(one ring matches)
        print(f"[{label}] {src}  length {n}  non-ACGT {n - tot}  "
              f"GC={pC + pG:.4f}  q={q:.5f}", flush=True)
        for k in KS:
            dt, (pos, cen, pairs) = t_of(expand, g, k)
            exp_c = (n - k + 1) * q ** (k // 2)
            rows.append(dict(genome=label, k=k, seconds=dt,
                             start_positions=len(pos), centres=cen,
                             pairs=pairs, expected_centres=exp_c,
                             obs_over_exp=cen / exp_c))
            print("   ", rows[-1], flush=True)
    save("palindromes_genome", rows)

    # ---- side-by-side summary, only when a real genome was analysed
    by = {(r["genome"], r["k"]): r for r in rows}
    if ("real", KS[0]) in by:
        print("\nobs/exp comparison (deficit = 1 - real/synthetic):")
        for k in KS:
            r, s = by[("real", k)], by[("synthetic", k)]
            print(f"  k={k:2d}  real {r['obs_over_exp']:.4f}   "
                  f"synthetic {s['obs_over_exp']:.4f}   "
                  f"deficit {1 - r['obs_over_exp'] / s['obs_over_exp']:+.3f}")

    # ---- (b) naive vs expansion on prefixes of the first genome
    g0 = genomes[0][1]
    pre = []
    for k in KS:
        for m in PREFIX_SIZES:
            s = g0[:m]
            dn, pr = t_of(naive_pairs, s, k)
            de, ex = t_of(expand, s, k)
            assert summarize_pairs(pr) == ex, (k, m)
            pre.append(dict(k=k, n=m, t_naive=dn, t_expand=de, pairs=len(pr)))
            print(pre[-1], flush=True)
    save("palindromes_prefix", pre)

    # ---- (c) worst-case families (k=4)
    wc = []
    for fam, mk in (("A^(n/2)T^(n/2)", lambda n: "A" * (n // 2) + "T" * (n // 2)),
                    ("(AT)^(n/2)",     lambda n: "AT" * (n // 2))):
        for m in WORST_SIZES:
            s = mk(m)
            dn, pr = t_of(naive_pairs, s, 4)
            de, ex = t_of(expand, s, 4)
            assert summarize_pairs(pr) == ex, (fam, m)
            wc.append(dict(family=fam, n=m, t_naive=dn, t_expand=de, pairs=len(pr)))
            print(wc[-1], flush=True)
    save("palindromes_worst", wc)


if __name__ == "__main__":
    main()

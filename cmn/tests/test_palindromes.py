import random, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from palindromes import naive_pairs, expand_pairs, expand, summarize_pairs

_C = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(_C)[::-1]          # independent oracle
def is_rc_pal(w):
    # a window containing a non-ACGT symbol (e.g. N) is never a palindrome
    return set(w) <= set("ACGT") and w == revcomp(w)
def oracle_pairs(t, k):
    return [(i, l) for i in range(len(t)) for l in range(k, len(t) - i + 1)
            if is_rc_pal(t[i:i+l])]

def test_gaattc():
    assert revcomp("GAATTC") == "GAATTC"
    for f in (naive_pairs, expand_pairs):
        assert f("GAATTC", 6) == [(0, 6)]
        assert f("GAATTC", 4) == [(0, 6), (1, 4)]
        assert f("GAATTC", 2) == [(0, 6), (1, 4), (2, 2)]
    assert expand("GAATTC", 4)[0] == [0, 1]

def test_edge_cases():
    for f in (naive_pairs, expand_pairs):
        assert f("", 2) == []
        assert f("A", 1) == []
        assert f("AT", 3) == []                    # k > n
        assert f("AT", 2) == [(0, 2)]
        assert f("NNNN", 2) == []                  # N is not its own complement
        assert f("gaattc", 6) == [(0, 6)]          # lower case accepted
        assert f("AATT", 3) == f("AATT", 4) == [(0, 4)]   # odd k rounds up
    for bad in (0, -1):
        for f in (naive_pairs, expand_pairs, expand):
            try: f("AT", bad); assert False
            except ValueError: pass

def test_random_against_oracle():
    rng = random.Random(1)
    for _ in range(1500):
        n = rng.randint(0, 40)
        al = rng.choice(["ACGT", "AT", "ACGTN", "AC"])
        t = "".join(rng.choice(al) for _ in range(n))
        k = rng.randint(1, 9)
        want = oracle_pairs(t, k)
        assert naive_pairs(t, k) == want, (t, k)
        assert expand_pairs(t, k) == want, (t, k)
        assert expand(t, k) == summarize_pairs(want), (t, k)

def test_AT_repeat_all_even_windows_are_palindromes():
    n, k = 20, 4
    t = "AT" * (n // 2)
    want = sum(n - l + 1 for l in range(k, n + 1, 2))
    assert len(expand_pairs(t, k)) == want == len(naive_pairs(t, k))

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok  ", name)

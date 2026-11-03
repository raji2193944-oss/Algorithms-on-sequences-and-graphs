import random, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from common_substring import naive, hashing

def brute_exists(x, y, k):
    return any(x[i:i+k] in y for i in range(len(x) - k + 1))

def valid(x, y, k, r):
    i, j = r
    return x[i:i+k] == y[j:j+k]

def test_edge_cases():
    for f in (naive, hashing):
        assert f("", "", 1) is None
        assert f("ABC", "", 1) is None
        assert f("ABC", "ABC", 4) is None          # k > n
        assert f("ABC", "ABC", 3) == (0, 0)        # k == n == m
        assert f("ABC", "XBC", 3) is None
        assert f("ABCD", "XXCD", 2) == (2, 2)
        assert f("A", "A", 1) == (0, 0)

def test_random_small_alphabets():
    rng = random.Random(0)
    for _ in range(3000):
        a = rng.choice([2, 3, 4, 20])
        al = "ACDEFGHIKLMNPQRSTVWY"[:a]
        n, m = rng.randint(0, 25), rng.randint(0, 25)
        x = "".join(rng.choice(al) for _ in range(n))
        y = "".join(rng.choice(al) for _ in range(m))
        k = rng.randint(1, 8)
        truth = brute_exists(x, y, k)
        for f in (naive, hashing):
            r = f(x, y, k)
            assert (r is not None) == truth, (f.__name__, x, y, k)
            if r is not None:
                assert valid(x, y, k, r)

def test_invalid_k():
    for f in (naive, hashing):
        for bad in (0, -3):
            try:
                f("AB", "AB", bad); assert False
            except ValueError:
                pass

def test_adversarial_instance_has_no_match():
    k = 7
    x = "A" * 60
    y = ("A" * (k - 1) + "C") * 9
    assert naive(x, y, k) is None and hashing(x, y, k) is None


if __name__ == "__main__":      # runs without pytest: python tests/test_common_substring.py
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok  ", name)

"""Task 2.2 -- reverse-complement palindromes.

A window s is a reverse-complement palindrome iff s == revcomp(s), i.e. for
every a:  s[a] == complement(s[len-1-a]).  Characters outside ACGT never match.
Note: such windows always have EVEN length (a base is never its own complement).

Output convention used throughout ("pairs"): the set of (i, l) with l >= k such
that t[i:i+l] is a reverse-complement palindrome.  "positions" = sorted {i}.
"""

COMP = bytearray(256)                       # 0 for anything that is not ACGT
for _a, _b in ("AT", "TA", "CG", "GC"):
    COMP[ord(_a)] = ord(_b)


def _prep(t):
    if isinstance(t, str):
        t = t.upper().encode("ascii", "replace")
    return t


def _check_k(k):
    if k < 1:
        raise ValueError("k must be >= 1")


def naive_pairs(t, k):
    """Every start i, every length l >= k, test the window outside-in.
    Worst case Theta(n^3); on random DNA Theta(n^2) (each test fails at once)."""
    _check_k(k)
    s = _prep(t)
    n = len(s)
    out = []
    for i in range(n):
        for l in range(k, n - i + 1):
            a, b = i, i + l - 1
            while a <= b and s[a] == COMP[s[b]]:
                a += 1
                b -= 1
            if a > b:                       # all rings matched (never for odd l)
                out.append((i, l))
    return out


def _radii(s):
    """Yield (c, r): centre c lies between s[c-1] and s[c]; r = maximal radius."""
    n = len(s)
    for c in range(1, n):
        r = 0
        lo, hi = c - 1, c
        while lo >= 0 and hi < n and s[lo] == COMP[s[hi]]:
            r += 1
            lo -= 1
            hi += 1
        yield c, r


def expand_pairs(t, k):
    """Centre expansion, listing every (i, l) pair.  Time O(n*k + P)."""
    _check_k(k)
    s = _prep(t)
    rho0 = (k + 1) // 2                     # smallest radius with length 2*rho >= k
    out = []
    for c, r in _radii(s):
        for rho in range(rho0, r + 1):
            out.append((c - rho, 2 * rho))
    out.sort()
    return out


def expand(t, k):
    """Centre expansion, positions only.  Returns (positions, n_centres, n_pairs).

    n_centres = number of centres whose maximal palindrome has length >= k
    n_pairs   = P, the number of (i, l) pairs.
    Start positions of a centre with radius r form the interval [c-r, c-rho0];
    a difference array merges the intervals in O(1) each, so the time is
    O(n + sum_c r_c)  (<= O(n*k + P))  and the memory is O(n)."""
    _check_k(k)
    s = _prep(t)
    n = len(s)
    rho0 = (k + 1) // 2
    diff = [0] * (n + 2)
    centres = pairs = 0
    for c, r in _radii(s):
        if r >= rho0:
            centres += 1
            pairs += r - rho0 + 1
            diff[c - r] += 1
            diff[c - rho0 + 1] -= 1
    positions, run = [], 0
    for i in range(n):
        run += diff[i]
        if run > 0:
            positions.append(i)
    return positions, centres, pairs


def summarize_pairs(pairs):
    """(positions, n_centres, n_pairs) computed from an explicit pair list."""
    pos = sorted({i for i, _ in pairs})
    centres = len({i + l // 2 for i, l in pairs})
    return pos, centres, len(pairs)

if __name__ == "__main__":
    t = "GAATTC"                    # EcoRI site, a length-6 RC palindrome
    print(naive_pairs(t, 6))        # [(0, 6)]
    print(expand_pairs(t, 6))       # [(0, 6)]
    print(expand(t, 6))             # ([0], 1, 1)

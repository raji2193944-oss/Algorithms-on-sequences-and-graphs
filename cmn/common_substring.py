"""Task 2.1 -- common substring of length exactly k.

Both functions return (i, j) with x[i:i+k] == y[j:j+k], or None.
"""


def naive(x, y, k):
    """For every window of x, compare against every window of y.

    The inner while-loop is written out character by character on purpose:
    it makes the Theta(k) cost of one window comparison visible in the code,
    which is what the O(n*m*k) proof counts.
    """
    if k <= 0:
        raise ValueError("k must be positive")
    n, m = len(x), len(y)
    for i in range(n - k + 1):
        for j in range(m - k + 1):
            l = 0
            while l < k and x[i + l] == y[j + l]:
                l += 1
            if l == k:
                return (i, j)
    return None


def hashing(x, y, k):
    """Store all k-mers of x in a dict (k-mer -> first position), scan y."""
    if k <= 0:
        raise ValueError("k must be positive")
    n, m = len(x), len(y)
    seen = {}
    for i in range(n - k + 1):
        seen.setdefault(x[i:i + k], i)      # O(k) to slice + O(k) to hash
    for j in range(m - k + 1):
        i = seen.get(y[j:j + k])            # O(k) to slice + O(k) to hash
        if i is not None:
            return (i, j)
    return None
if __name__ == "__main__":
    print(naive("ACDEFG", "GGDETT", 2))
    print(hashing("ACDEFG", "GGDETT", 2))

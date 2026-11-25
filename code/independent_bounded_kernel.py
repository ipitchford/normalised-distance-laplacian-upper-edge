"""Check all complete graphs of orders 4, 5 and 6 with edge weights 1 or 3."""
from fractions import Fraction
from itertools import combinations
from math import floor, ceil

import numpy as np



def require(condition):
    """Fail loudly even under python -O (no reliance on assert)."""
    if not condition:
        raise SystemExit("CHECK FAILED")

def balanced_value(n):
    def f(a):
        return Fraction(a - 1, (a - 1) + 3 * (n - a))
    return f(floor(n / 2)) + f(ceil(n / 2))


def check(n):
    edges = list(combinations(range(n), 2))
    cases = 2 ** len(edges)
    best = float("inf")
    for first in range(0, cases, 256):
        masks = np.arange(first, min(first + 256, cases))
        weights = np.zeros((len(masks), n, n))
        for index, (i, j) in enumerate(edges):
            weight = 1 + 2 * ((masks >> index) % 2)
            weights[:, i, j] = weight
            weights[:, j, i] = weight
        t = weights.sum(axis=2)
        normalised = weights / np.sqrt(t[:, :, None] * t[:, None, :])
        smallest = np.linalg.eigvalsh(normalised)[:, 0]
        best = min(best, 1 + smallest.min())
    expected = balanced_value(n)
    require(abs(float(expected) - best) < 1e-12)
    print(n, cases, expected, f"{best:.12f}")


if __name__ == "__main__":
    print("order cases exact_block_bound smallest_computed_deficit")
    for n in (4, 5, 6):
        check(n)

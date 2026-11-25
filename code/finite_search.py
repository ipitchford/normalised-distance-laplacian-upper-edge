"""Finite evidence for the normalised distance Laplacian upper edge (distance power one).

Three computations, all reported in Section 11 (Finite computations) of the manuscript:

1. exhaustive  -- every connected graph of order 4..9 (nauty geng -c) and every tree of
   order 6..TREE_MAX (nauty gentreeg). For each order the script finds the minimum deficit
   eps(G) = 1 + lambda_min(T^{-1} D) and compares it with the best member of the
   clique-path-clique family (all graphs) or of the double-broom family (trees: two stars
   whose centres are joined by a path).
2. optimum     -- the best member found (ternary search in the path length, then a local
   window) among balanced and near-balanced members of these two
   families for n up to 10^6, computed on the equitable quotient (the vertices of each
   clique or star, other than its attachment vertex, are interchangeable).
3. self-tests  -- the quotient eigenvalue is compared with a full distance-matrix
   eigenvalue for small n.

Floating point (LAPACK) is used for eigenvalues. Declared tolerance: two deficits are
treated as equal when they differ by less than 1e-10; the printed table values are
compared across platforms to 1e-9 relative. Nothing here is an ingredient of a proof.

Usage: python3 finite_search.py [--tree-max 20] [--output ../results/finite_search.json]
Requires numpy, scipy and the nauty programs geng and gentreeg on PATH (or NAUTY_BIN).
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

TOL = 1e-10


def require(condition, message="CHECK FAILED"):
    if not condition:
        raise SystemExit(message)


def nauty(prog):
    base = os.environ.get("NAUTY_BIN")
    path = os.path.join(base, prog) if base else shutil.which(prog)
    require(path is not None and os.path.exists(path), f"nauty program {prog} not found")
    return path


def decode_graph6(line):
    """Decode graph6 (no header) to an adjacency matrix (n <= 62)."""
    data = [c - 63 for c in line.encode()]
    n = data[0]
    bits = []
    for v in data[1:]:
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    a = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                a[i, j] = a[j, i] = 1
            idx += 1
    return a


def decode_parent_array(line):
    """Decode a gentreeg -p line (1-indexed parents, root parent 0)."""
    parents = [int(x) for x in line.split()]
    n = len(parents)
    a = np.zeros((n, n), dtype=np.int8)
    for child, parent in enumerate(parents):
        if parent:
            a[child, parent - 1] = a[parent - 1, child] = 1
    require(int(a.sum()) == 2 * (n - 1), "bad parent array")
    return a


def deficit_from_adjacency(a):
    d = shortest_path(a, unweighted=True, directed=False)
    t = d.sum(axis=1)
    return 1.0 + np.linalg.eigvalsh(d / np.sqrt(np.outer(t, t)))[0]


def kpk_adjacency(n1, inner, n3):
    """Cliques K_n1, K_n3 joined by a path with `inner` internal vertices."""
    n = n1 + inner + n3
    a = np.zeros((n, n), dtype=np.int8)
    a[:n1, :n1] = 1
    a[n1:n1 + n3, n1:n1 + n3] = 1
    np.fill_diagonal(a, 0)
    chain = [0] + list(range(n1 + n3, n)) + [n1]
    for u, v in zip(chain, chain[1:]):
        a[u, v] = a[v, u] = 1
    return a


def broom_adjacency(leaves1, inner, leaves2):
    """Stars with leaves1, leaves2 leaves whose centres are joined by `inner` path vertices."""
    n = leaves1 + leaves2 + inner + 2
    a = np.zeros((n, n), dtype=np.int8)
    chain = list(range(inner + 2))
    for u, v in zip(chain, chain[1:]):
        a[u, v] = a[v, u] = 1
    nxt = inner + 2
    for centre, count in ((0, leaves1), (inner + 1, leaves2)):
        for _ in range(count):
            a[centre, nxt] = a[nxt, centre] = 1
            nxt += 1
    return a


def family_best_small(n, tree):
    best = (math.inf, None)
    if tree:
        for inner in range(0, n - 1):
            for l1 in range(0, n - 1 - inner):
                l2 = n - 2 - inner - l1
                if l2 < 0:
                    continue
                e = deficit_from_adjacency(broom_adjacency(l1, inner, l2))
                if e < best[0]:
                    best = (e, (l1, inner, l2))
    else:
        for inner in range(0, n - 1):
            for n1 in range(1, n - inner):
                n3 = n - inner - n1
                if n3 < 1:
                    continue
                e = deficit_from_adjacency(kpk_adjacency(n1, inner, n3))
                if e < best[0]:
                    best = (e, (n1, inner, n3))
    return best


def exhaustive(tree_max):
    rows = []
    geng, gentreeg = nauty("geng"), nauty("gentreeg")
    for n in range(4, 10):
        out = subprocess.run([geng, "-c", "-q", str(n)], capture_output=True, text=True, check=True)
        graphs = out.stdout.split()
        values = sorted(deficit_from_adjacency(decode_graph6(g)) for g in graphs)
        fam, params = family_best_small(n, tree=False)
        rows.append({"class": "connected graphs", "n": n, "count": len(graphs),
                     "min_deficit": values[0], "family_best": fam, "family_params": params,
                     "second_smallest": values[1], "match": abs(values[0] - fam) < TOL})
        print(f"graphs n={n} count={len(graphs)} min={values[0]:.10f} KPK={fam:.10f} {params}", flush=True)
    for n in range(6, tree_max + 1):
        out = subprocess.run([gentreeg, "-p", "-q", str(n)], capture_output=True, text=True, check=True)
        trees = out.stdout.strip().split("\n")
        best = math.inf
        second = math.inf
        for g in trees:
            e = deficit_from_adjacency(decode_parent_array(g))
            if e < best:
                best, second = e, best
            elif e < second:
                second = e
        fam, params = family_best_small(n, tree=True)
        rows.append({"class": "trees", "n": n, "count": len(trees), "min_deficit": best,
                     "family_best": fam, "family_params": params, "second_smallest": second,
                     "match": abs(best - fam) < TOL})
        print(f"trees  n={n} count={len(trees)} min={best:.10f} broom={fam:.10f} {params}", flush=True)
    return rows


def quotient_deficit(n1, edges, n3, tree):
    """Exact deficit of a balanced-type member via its equitable quotient.

    Vertex classes: the n1-1 non-attachment block vertices, the edges+1 path vertices
    (attachment vertices included), and the n3-1 other block vertices. `tree` selects
    stars (internal block distance 2) instead of cliques (internal distance 1).
    """
    din = 2.0 if tree else 1.0
    r = edges
    sizes = np.array([n1 - 1] + [1] * (r + 1) + [n3 - 1], dtype=float)
    m = len(sizes)
    pos = np.arange(r + 1, dtype=float)
    dist = np.zeros((m, m))
    dist[1:-1, 1:-1] = np.abs(pos[:, None] - pos[None, :])
    dist[0, 1:-1] = dist[1:-1, 0] = 1 + pos
    dist[-1, 1:-1] = dist[1:-1, -1] = 1 + (r - pos)
    dist[0, -1] = dist[-1, 0] = r + 2
    dist[0, 0] = dist[-1, -1] = din
    pair = np.outer(sizes, sizes) * dist
    pair[0, 0] = sizes[0] * (sizes[0] - 1) * din
    pair[-1, -1] = sizes[-1] * (sizes[-1] - 1) * din
    keep = sizes > 0
    pair, sizes = pair[np.ix_(keep, keep)], sizes[keep]
    t = pair.sum(axis=1) / sizes
    s = pair / np.sqrt(np.outer(t * sizes, t * sizes))
    return 1.0 + eigh(s, eigvals_only=True, subset_by_index=[0, 0])[0]


def family_optimum(n, tree):
    """Minimise over the number of path edges r and splits within 6 of balance."""
    def balanced(r):
        rest = n - (r - 1)
        n1 = rest // 2
        return quotient_deficit(n1, r, rest - n1, tree)
    scale = math.sqrt((6 if tree else 3) * n)
    lo, hi = max(1, int(0.6 * scale)), int(1.5 * scale) + 3
    while hi - lo > 3:
        m1, m2 = lo + (hi - lo) // 3, hi - (hi - lo) // 3
        if balanced(m1) < balanced(m2):
            hi = m2
        else:
            lo = m1
    r_best = min(range(lo, hi + 1), key=balanced)
    best = (balanced(r_best), r_best, (n - (r_best - 1)) // 2)
    for dr in range(-2, 3):
        r = r_best + dr
        rest = n - (r - 1)
        for shift in range(-6, 7):
            n1 = rest // 2 + shift
            if r < 1 or n1 < 2 or rest - n1 < 2:
                continue
            e = quotient_deficit(n1, r, rest - n1, tree)
            if e < best[0] - 1e-15:
                best = (e, r, n1)
    return best


def optimum_table():
    rows = []
    for tree in (False, True):
        b = 2.0 if tree else 1.0
        limit = 4 * math.sqrt(b / 3)
        for n in (100, 1000, 10**4, 10**5, 10**6):
            e, r, n1 = family_optimum(n, tree)
            rows.append({"family": "double broom" if tree else "clique-path-clique", "n": n,
                         "path_edges": r, "block_sizes": [n1, n - (r - 1) - n1],
                         "deficit": e, "sqrt_n_deficit": math.sqrt(n) * e,
                         "r_over_sqrt_n": r / math.sqrt(n),
                         "second_order_coefficient": (math.sqrt(n) * e - limit) * math.sqrt(n),
                         "limit": limit})
            print(f"{rows[-1]['family']:19s} n={n:>8d} r={r:5d} sqrt(n)eps={math.sqrt(n)*e:.6f} "
                  f"r/sqrt(n)={r/math.sqrt(n):.4f} n*(eps-K/sqrt n)={rows[-1]['second_order_coefficient']:.4f}",
                  flush=True)
    return rows


def self_test():
    for tree in (False, True):
        for n1, r, n3 in ((3, 2, 4), (5, 4, 5), (7, 6, 6), (9, 3, 8)):
            q = quotient_deficit(n1, r, n3, tree)
            if tree:
                full = deficit_from_adjacency(broom_adjacency(n1 - 1, r - 1, n3 - 1))
            else:
                full = deficit_from_adjacency(kpk_adjacency(n1, r - 1, n3))
            require(abs(q - full) < TOL, f"quotient self-test failed {tree} {(n1, r, n3)}")
    print("PASS quotient self-test (8 cases)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree-max", type=int, default=20)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    self_test()
    ex = exhaustive(args.tree_max)
    require(all(row["match"] for row in ex), "a family member is not extremal")
    print("PASS exhaustive: family members attain the minimum deficit at every order checked", flush=True)
    opt = optimum_table()
    result = {"tolerance": TOL, "python": sys.version.split()[0], "numpy": np.__version__,
              "exhaustive": ex, "family_optimum": opt}
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(result, fh, indent=1, default=lambda o: o.item())
    print("PASS")


if __name__ == "__main__":
    main()

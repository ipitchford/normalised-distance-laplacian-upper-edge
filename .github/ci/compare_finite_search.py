# SPDX-License-Identifier: MIT
"""Compare a fresh results/finite_search.json with the recorded one under the declared policy.

usage: compare_finite_search.py NEW RECORDED [--subset]

Policy (declared before tagging; the same in README.md, ENVIRONMENT.txt and REPLAY_RECEIPT.md):
- the top-level fields "python" and "numpy" (versions) are ignored;
- floating-point values (LAPACK eigenvalues and quantities derived from them) agree within 1e-9 relative, except
  "second_order_coefficient" (an eigenvalue difference multiplied by n), compared within 1e-4 absolute;
- every key, integer, string, boolean (including every "match") and list length is identical.
--subset: NEW may come from a shorter run (smaller --tree-max); each of its exhaustive rows is compared with the
recorded row of the same class and order, and every other part must match in full.
Exits 1 on any difference.
"""
import json
import sys


def walk(a, b, path, bad, counts):
    if isinstance(a, bool) or isinstance(b, bool) or not (isinstance(a, float) or isinstance(b, float)):
        if type(a) is not type(b):
            bad.append(f'{path}: type {type(a).__name__} vs {type(b).__name__}')
        elif isinstance(b, dict):
            if list(a) != list(b):
                bad.append(f'{path}: keys {list(a)} vs {list(b)}')
            else:
                for k in b:
                    walk(a[k], b[k], f'{path}.{k}', bad, counts)
        elif isinstance(b, list):
            if len(a) != len(b):
                bad.append(f'{path}: length {len(a)} vs {len(b)}')
            else:
                for i, (x, y) in enumerate(zip(a, b)):
                    walk(x, y, f'{path}[{i}]', bad, counts)
        elif a != b:
            bad.append(f'{path}: {a!r} vs {b!r}')
        else:
            counts['exact'] += 1
        return
    if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        bad.append(f'{path}: type {type(a).__name__} vs {type(b).__name__}')
    elif path.endswith('.second_order_coefficient'):
        if abs(a - b) > 1e-4:
            bad.append(f'{path}: {a!r} vs {b!r} beyond 1e-4 absolute')
        else:
            counts['float'] += 1
    elif abs(a - b) > 1e-9 * max(abs(a), abs(b)):
        bad.append(f'{path}: {a!r} vs {b!r} beyond 1e-9 relative')
    else:
        counts['float'] += 1


def main():
    new = json.load(open(sys.argv[1], encoding='utf-8'))
    old = json.load(open(sys.argv[2], encoding='utf-8'))
    subset = '--subset' in sys.argv[3:]
    bad, counts = [], {'exact': 0, 'float': 0}
    for d in (new, old):
        for k in ('python', 'numpy'):
            d.pop(k, None)
    if subset:
        rec = {(r['class'], r['n']): r for r in old['exhaustive']}
        rows = new.pop('exhaustive')
        old.pop('exhaustive')
        if not rows:
            bad.append('no exhaustive rows')
        for r in rows:
            key = (r.get('class'), r.get('n'))
            if key not in rec:
                bad.append(f'no recorded row for {key}')
            else:
                walk(r, rec[key], f'$.exhaustive{list(key)}', bad, counts)
    walk(new, old, '$', bad, counts)
    if bad:
        print('DIFFERENCES:')
        for m in bad:
            print(' ', m)
        sys.exit(1)
    print(f'OK: {counts["exact"]} values identical, {counts["float"]} floats within tolerance; "python" and "numpy" ignored'
          + (' (subset of exhaustive rows)' if subset else ''))


if __name__ == '__main__':
    main()

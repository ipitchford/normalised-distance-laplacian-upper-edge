# Replay receipt (v0.1.0-candidate, 25 November 2026)

Release replay by the publisher from fresh copies of the approved candidate (code and requirements only; results and
figure regenerated), working directory = `code/`, macOS arm64, CPython 3.12.11 with numpy 2.3.5, scipy 1.17.0 and
sympy 1.14.0 (printed inside the runs; see [ENVIRONMENT.txt](ENVIRONMENT.txt)).

| Run | Command (from `code/`) | Exit | Output vs recorded |
|---|---|---|---|
| verify -O | `uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with sympy==1.14.0 python3 -I -B -O verify.py` | 0 | identical |
| independent_trial_quotient -O | same environment | 0 | identical |
| independent_bounded_kernel -O | same environment | 0 | identical |
| finite_search -O | `uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 python3 -I -B -O finite_search.py --tree-max 20 --output ../results/finite_search.json` (1 CPU, 505 s) | 0 | finite_search.txt and finite_search.json identical byte for byte; three PASS lines |
| figure | `uv run --python 3.12 --with numpy==2.3.5 --with matplotlib==3.10.8 python3 -I -B figure.py` | 0 | PDF bytes differ only in the creation date; rendering identical at 100 dpi |
| Negative control | verify.py with the check `require(poly.subs(x,1)==0)` changed to `==1`, under -O | 1 | `CHECK FAILED` |

stderr was empty for every passing run. The lead's earlier clean replay (24 November 2026, CPython 3.12.11) also
produced byte-identical outputs. After the replay only the title-page DOI line of manuscript.tex, the rebuilt PDF and
the release records changed; `code/`, `results/` and `figures/` are byte-identical to the replayed candidate.

Comparison policy for other platforms (declared before tagging; also in CI, README.md and ENVIRONMENT.txt): the three
quick text outputs byte-identical; finite_search.json compared as JSON with versions ignored, floats within 1e-9
relative and `second_order_coefficient` within 1e-4 absolute, everything else identical.

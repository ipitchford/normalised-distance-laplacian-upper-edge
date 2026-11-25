# The sharp upper edge of the normalised distance Laplacian

**How close the distance Laplacian can get to two.** Unrefereed research candidate, version 0.1.0-candidate.
DOI [10.5281/zenodo.24540994](https://doi.org/10.5281/zenodo.24540994) (all versions:
[10.5281/zenodo.24540993](https://doi.org/10.5281/zenodo.24540993)). Release page:
https://evidencepress.org/releases/normalised-distance-laplacian-upper-edge/

## What this is about

**The picture.** Take any connected network and a walker that jumps from node to node, preferring nodes that are far
away: the chance of jumping to a node is proportional to its distance. In some networks this walker ends up bouncing
back and forth between two distant regions, almost like a pendulum. One number measures how close it comes to perfect
back-and-forth: the largest eigenvalue of the *normalised distance Laplacian*. It always lies between 1 and 2, can come
close to 2, and never reaches it. The question is: over all networks with n nodes, how close to 2 can it get, and which
networks get closest? The gap from 2 is called the deficit.

**What was known.** Reinhart introduced the matrix (Special Matrices, 2021), computed the best networks with up to ten
nodes, and conjectured that the winners are two complete clusters ("cliques") joined by a path, with the value tending
to 2. Johnston and Tait (Linear Algebra and its Applications, 2024) proved that the value tends to 2. Their results put
the smallest deficit between a constant times 1/n and about 2.83/√n (the second from a two-cliques-and-a-path
construction, constant 2√2), proved the order 1/√n under extra conditions on the eigenvector, and asked for the true
order of the gap.

**What is new.**

- The smallest deficit over networks with n nodes is 4/√(3n) to first order, about 2.31/√n, with no extra assumptions.
- Every network that gets this close looks, at large scale, like a path of length about √(3n) with half of the nodes
  packed into a nearly complete cluster at each end. The eigenvector that measures the back-and-forth falls linearly
  along the path; that linear profile is what beats the earlier 2.83 construction, where it was zero on the path.
- The same method gives sharp constants when distance is raised to any fixed positive power, for trees (4√(2/3), about
  3.27, at power one) and for networks that avoid a fixed sub-pattern containing a cycle. It also describes the limiting
  eigenvector and the rest of the spectrum, and treats slowly growing and bounded distance weights.
- New computer evidence: among all trees with 6 to 20 nodes, the best tree is always two stars whose centres are joined
  by a path. The paper states this for every n as a conjecture.

Figure 1 of the paper shows the balance behind the constant: with a path of length s√n, the cost is about 2/s from
pairs inside the two dense ends and 2s/3 from pairs between the path and the ends; the total 2/s + 2s/3 is smallest at
s = √3, where it equals 4/√3.

**Idea of the proof.** Look at the eigenvector at its highest and lowest nodes. Because the walk almost flips the sign
of the eigenvector, the jump probabilities from these two nodes can share very little; this is the classical
Doeblin–Dobrushin bound, applied to these two rows (no novelty is claimed for the bound itself). It forces almost all
nodes to crowd near one of the two. Each node on a shortest path between the two crowds then pays a definite price,
whatever eigenvector value it takes, and balancing that price against the cost inside the crowds gives the constant.

**Limits.** The result is first order. Networks that differ from two cliques joined by a path by a small fraction of
edges score the same to this order, so Reinhart's exact conjecture stays open. The next term is observed only
numerically (about 4.26/n at a million nodes for the best clique–path–clique network, still drifting). The tree
statement is a conjecture backed by computation up to 20 nodes. This is an unrefereed candidate; its assurance is
internal AI-assisted review and computational checks only (see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)).

**Credit.** The question, the barbell construction and the conditional order are due to Johnston and Tait; the matrix,
the conjecture and the small-order table to Reinhart.

## Where to read next

- The paper: [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex), figure
  [figures/distance-balance.pdf](figures/distance-balance.pdf)). Section 1 is a plain introduction with what was known,
  what is new, what is not settled and a roadmap.
- Machine-readable map of claims and evidence: [AI_INDEX.md](AI_INDEX.md).
- Prior art, sources read and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json).
- Review summary (internal, AI-assisted): [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).

## Contents

| Path | What it is |
|---|---|
| `manuscript.tex`, `manuscript.pdf`, `figures/distance-balance.pdf` | The paper (16 pages) and its figure |
| `code/verify.py` | Exact characteristic polynomials of small paths; exact integer trial quotients of the constructions up to n = 10^8; full distance-matrix spectra for powers 0.5, 1, 2, 3 at n = 211; logarithmic and d log(1+d) weights; the bounded-weight inequality on 180 weight matrices |
| `code/finite_search.py` | All connected graphs of order 4 to 9 (nauty `geng -c`) and all trees of order 6 to 20 (`gentreeg`); best family members up to n = 10^6 on the equitable quotient; quotient self-test |
| `code/independent_trial_quotient.py`, `code/independent_bounded_kernel.py` | Separate recomputations of the trial quotients and of the bounded-weight inequality (all weight patterns in {1, 3} for n = 4, 5, 6) |
| `code/figure.py` | Draws Figure 1 |
| `results/` | Recorded outputs of the four verification scripts |
| `refs.bib` | Bibliography (Crossref and arXiv records) |
| `ENVIRONMENT.txt`, `REPLAY_RECEIPT.md`, `.github/` | Replay environment, receipt and CI |
| `MANIFEST.sha256` | SHA-256 of every other file |

## Replay

Requirements: Python 3.12, the pinned packages in `requirements.txt`, and for `finite_search.py` the nauty programs
`geng` and `gentreeg` on PATH (or in the folder named by `NAUTY_BIN`). Run from `code/` (scripts write `../results` and
`../figures`):

```sh
cd code
python3 -I -B -O verify.py > ../results/verify.txt
python3 -I -B -O independent_trial_quotient.py > ../results/independent_trial_quotient.txt
python3 -I -B -O independent_bounded_kernel.py > ../results/independent_bounded_kernel.txt
python3 -I -B -O finite_search.py --tree-max 20 --output ../results/finite_search.json > ../results/finite_search.txt
python3 -I -B figure.py
```

The first three take seconds; `finite_search.py --tree-max 20` takes about ten minutes on one core (`--tree-max 14`
takes a few minutes and is a quick smoke test). Every check raises an error and exits non-zero on failure, also under
`python -O`; `finite_search.py` prints `PASS` three times (self-test, exhaustive, final).

Comparison policy (declared before tagging; the same in CI, ENVIRONMENT.txt and REPLAY_RECEIPT.md):

- `results/verify.txt`, `results/independent_trial_quotient.txt` and `results/independent_bounded_kernel.txt`:
  identical byte for byte (they are deterministic; exact or rational where the claims depend on it).
- `results/finite_search.json`: compared as JSON with `.github/ci/compare_finite_search.py`. The fields `python` and
  `numpy` (versions) are ignored; floating-point values (LAPACK eigenvalues) agree within 1e-9 relative, except the
  derived field `second_order_coefficient`, which multiplies an eigenvalue difference by n and is compared within 1e-4
  absolute (the precision printed in the paper); every integer, string, boolean (including every `match`) and list is
  identical. `results/finite_search.txt` prints rounded values of the same numbers and must contain the three `PASS` lines.
- `figures/distance-balance.pdf` is regenerated with a new creation date; its rendering is identical.

GitHub CI runs the three quick scripts under `-O`, compares them byte for byte, runs a negative control, and checks the
manifest. The full `finite_search.py` run is kept off the blocking path; it was replayed on the release host (see
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)).

## Licence and citation

Text, the paper, the figure, recorded outputs and records: CC0-1.0. Code: MIT. See [LICENSE.md](LICENSE.md) and
[CITATION.cff](CITATION.cff). Creator: Anonymous.

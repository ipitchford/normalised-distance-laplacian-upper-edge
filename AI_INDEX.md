# AI index: The sharp upper edge of the normalised distance Laplacian (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction: the question, what was known, what is new, what is not settled, idea of the proof, roadmap).
DOI 10.5281/zenodo.24540994 (concept 10.5281/zenodo.24540993). Status: unrefereed candidate; internal AI-assisted
review and computational checks only.

## Normalisation

G connected simple graph, n vertices, d(u,v) graph distance, D distance matrix, T = diag(row sums) (transmissions).
P = T^-1 D (walk that jumps with probability proportional to distance). Normalised distance Laplacian (Reinhart)
L = I - T^-1/2 D T^-1/2; eigenvalues of L are 1 - lambda for eigenvalues lambda of P. rho(G) = largest eigenvalue of L
= 1 - lambda_min(P), in (1, 2). Deficit eps(G) = 2 - rho(G) = 1 + lambda_min(P). M_n = max rho over connected graphs of
order n. General weights: d^alpha (alpha > 0 fixed), h_n = n^(1/(1+alpha)); classes: all graphs (b = 1),
K_{kappa+1}-free (b = 1 + (2^alpha - 1)/kappa), trees (b = 2^alpha). C_alpha = 4 int_0^1 z^a(1-z)^a/(z^a+(1-z)^a) dz,
s_* = (2 alpha b / C_alpha)^(1/(1+alpha)), K_{alpha,b} = 2b s_*^-alpha + C_alpha s_* = min_s (2b s^-alpha + C_alpha s).
At alpha = 1: C_1 = 2/3, s_* = sqrt(3b), K_{1,b} = 4 sqrt(b/3).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A (main) | lim sqrt(n) (2 - M_n) = 4/sqrt3 (about 2.3094) | Theorem 1 | written proof (lower bound Prop. 5, construction Prop. 6); answers the order-of-magnitude question of Johnston-Tait (arXiv:2302.11459, end of Section 4) without their eigenvector conditions (their Theorem 4.7) |
| B | Any sequence with sqrt(n) eps(G_n) -> 4/sqrt3: diam = (sqrt3 + o(1)) sqrt n; two disjoint sets of (1/2 + o(1)) n vertices with induced edge density -> 1; measured Gromov-Hausdorff limit [0,1] with measure (delta_0 + delta_1)/2 | Theorem 1(a)-(c), Section 7 | written proof; first order only: allows o(n^2) edge changes in the ends and branches of length o(sqrt n) |
| C | General: lim h_n^alpha min eps = K_{alpha,b} for every fixed alpha > 0 and each of the three classes; near-extremal structure (shortest path of length (s_* + o(1)) h_n, two halves near its ends with mean d^alpha -> b) | Theorem 2, Sections 4-7 | written proof |
| D | Trees at alpha = 1: constant 4 sqrt(2/3) (about 3.2660); K_{kappa+1}-free: 4 sqrt((1+1/kappa)/3); graphs excluding a fixed subgraph F with a cycle | Theorem 2, Corollary 8 (Section 8) | written proof; Erdos-Stone and Turan via [es] |
| E | Limiting shape of the extremal eigenvector; all eigenvalues of L other than 0 and rho collapse to 1 in mean square | Corollary 9 (Section 9) | written proof |
| F | Regularly varying weights and bounded weights (contrast factor via Birkhoff's bound) | Theorem 10, Corollary 11 (Section 10) | written proof; Birkhoff bound credited via Gaubert-Qu Theorem 5.6 |
| G | Row-overlap inequality sum_x min(W_ix, W_jx) <= eps max(t_i, t_j) at the extremes i, j of an extremal eigenvector | Lemma 3 (Section 4) | written proof; this is the classical Doeblin-Dobrushin ergodicity-coefficient bound (Gaubert-Qu eqs (1)-(2)); novelty claimed only for applying it to these two rows |
| H | Computation: clique-path-clique graphs attain the minimum deficit among all connected graphs of order 4..9 (reproduces Reinhart Table 2.1 for these orders) | Section 11 | finite_search.py, exhaustive via nauty geng -c; floating point, tolerance 1e-10 |
| I | Computation: two stars whose centres are joined by a path attain the minimum deficit among all trees of order 6..20 (Table 2); paths for 6..9, both stars with >= 2 leaves from 10 | Section 11, Table 2; Conjecture 12 (all n) | finite_search.py, exhaustive via nauty gentreeg; conjecture not proved |
| J | Observation: best family members up to n = 10^6 (Table 3): sqrt(n) eps = 2.305138 (cliques), 3.257465 (stars) at 10^6; r/sqrt n -> sqrt3, sqrt6; second-order coefficient about 4.26 and 8.52 at 10^6, still increasing | Section 11, Table 3; open problem 2 | finite_search.py (ternary search plus window; values are upper bounds, deficits of explicit graphs); no claim beyond the observation |

## Evidence

| Check | Script | Output | Arithmetic |
|---|---|---|---|
| Exact characteristic polynomials on small paths; exact integer trial quotients E/H for both constructions at n = 10^2, 10^4, 10^6, 10^8 (2.308974489 and 3.265133160 at 10^8; limits 2.309401077, 3.265986324); direct integer energy checks for 260 parameter triples; full spectra at n = 211 for alpha in {0.5, 1, 2, 3}; log(1+d) and d log(1+d) weights; bounded-weight inequality on 180 matrices | [code/verify.py](code/verify.py) | [results/verify.txt](results/verify.txt) | exact (SymPy, integers, fractions) where claims depend on it; NumPy/SciPy floats for spectra |
| Trial quotients recomputed by separate code | [code/independent_trial_quotient.py](code/independent_trial_quotient.py) | [results/independent_trial_quotient.txt](results/independent_trial_quotient.txt) | exact |
| Bounded-weight inequality: all 2^C(n,2) weight patterns in {1,3}, n = 4, 5, 6 (64, 1024, 32768 cases; bounds 2/7, 7/20, 4/11 attained) | [code/independent_bounded_kernel.py](code/independent_bounded_kernel.py) | [results/independent_bounded_kernel.txt](results/independent_bounded_kernel.txt) | exact bound, float eigenvalues |
| H, I, J above; quotient self-test (8 cases) | [code/finite_search.py](code/finite_search.py) | [results/finite_search.txt](results/finite_search.txt), [results/finite_search.json](results/finite_search.json) | LAPACK floats; tolerance 1e-10 |
| Figure 1 (construction and the cost curve 2/s + 2s/3, minimum 4/sqrt3 at s = sqrt3) | [code/figure.py](code/figure.py) | [figures/distance-balance.pdf](figures/distance-balance.pdf) | closed forms |

Numbers quoted on the release page and in its audio: 4/sqrt3 about 2.31, sqrt3 about 1.73, 4 sqrt(2/3) about 3.27,
2 sqrt2 about 2.83 (earlier construction), trees with 6 to 20 vertices, graphs with 4 to 9 vertices, about 4.26 at a
million vertices. All are proved in the paper or printed by the scripts.

## Limits

First order only: the theorem does not decide Reinhart's conjecture that the maximiser at each order is a
clique-path-clique graph. The second-order term is observed numerically only. The tree statement (Conjecture 12) is
supported by computation for 6 to 20 vertices only. Power weights with fixed alpha; uniformity in alpha is open.
Open problems: Section 12.

## Replay and records

Replay commands and the comparison policy (three text outputs byte-identical; finite_search.json compared as JSON with
versions ignored, floats within 1e-9 relative and second_order_coefficient within 1e-4 absolute):
[README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml) with
[.github/ci/compare_finite_search.py](.github/ci/compare_finite_search.py). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
Prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib).
Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).

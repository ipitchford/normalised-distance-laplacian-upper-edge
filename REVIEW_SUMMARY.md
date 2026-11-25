# Review summary (internal, AI-assisted)

This is a rating-free public summary of the internal checking of version 0.1.0-candidate. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (24 November 2026)

1. **Combined five-role review panel** (one Claude subagent taking the roles of editor, methodology and proof rigour,
   spectral and extremal graph theory, a reader from outside the field, and devil's advocate) on the first draft, with
   the cited sources. The methodology role re-derived every step of the lower bound, the construction, the stability
   argument, the spectrum corollary and the weight results; all held. The panel recomputed the constants and ran its own
   exhaustive and large-order computations. It found no mathematical error and no critical issue. Its major items
   concerned exposition, attribution and the stated scope, not correctness.
2. **One revision round with extensions.** A new context-first introduction (the question, what was known, what is new,
   what is not settled, the idea of the proof, a roadmap); a notation table and an opening sentence for every section;
   the core of the lower bound rewritten as a proposition with four numbered claims and explicit error terms; a map from
   each clause of the main theorems to its proof; attribution of the row-overlap inequality to the Doeblin–Dobrushin
   ergodicity coefficient and of the bounded-weight factor to Birkhoff's bound (both via Gaubert–Qu); the scope
   statement in the abstract, introduction and discussion. Extensions: exhaustive computations over all connected
   graphs of order 4 to 9 and all trees of order 6 to 20; the tree conjecture (Conjecture 12); numerical second-order
   data up to a million vertices, recorded as an observation only. Every `assert` in the scripts was replaced by an
   explicit check that exits non-zero, also under `python -O`.
3. **Verification re-review and line-by-line proof audit** of the rewritten text (Sections 3 to 12). No critical or
   major issue. The minor items were applied: script outputs regenerated into `results/`; "best member found" in place
   of "exact optimum" for the ternary search; the Kumar–Patra bound stated correctly; a one-line justification of the
   eigenvalue form of the Dobrushin bound; notation clashes removed; Turán credited via the introduction of Erdős–Stone.
   One audit item was not applied, with evidence: the auditor thought Johnston and Tait give no constant for their
   construction, but their arXiv version states it explicitly (c = 2^{3/2} + o(1), after Proposition 4.2), and the text
   now quotes that remark.
4. **Clean replay** of all scripts from fresh copies, by the lead and again by the publisher, with a negative control
   (see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)).

No Astra/Codex second-model review took place: that route reported the model as unsupported for the account in this
batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| "The key inequality is classical." | Agreed. The row-overlap bound is the Doeblin–Dobrushin ergodicity-coefficient bound, credited via Gaubert–Qu (equations (1)–(2)). The paper claims no novelty for it; what is new is applying it to the two rows at the extremes of the eigenvector. |
| "First-order stability says nothing about Reinhart's exact conjecture." | Agreed and stated in the abstract, Section 1 ("What is not settled") and Section 12: graphs that differ from two cliques joined by a path by o(n^2) edges inside the ends, or by short branches, have the same leading deficit. |
| "Is the improvement over Johnston–Tait real?" | Their construction has constant 2√2 (stated in their arXiv version after Proposition 4.2); the linear eigenvector profile on the path gives 4/√3, and the lower bound shows 4/√3 is sharp. |
| "Hard to read." | Restructured as above, with a notation table and a proof map. |
| "Is the computational evidence used in the proofs?" | No. Section 11 says so; the computations confirm the constants, reproduce Reinhart's table for 4 to 9 vertices, and support the tree conjecture. Large-order values come from a ternary search and are upper bounds (deficits of explicit graphs). |
| Not addressed | Dobrushin's and Seneta's original texts, Bingham–Goldie–Teugels and a measured Gromov–Hausdorff reference are not cited directly; the paper cites only what was read, and credits the bounds through Gaubert–Qu. Two works cited in the first draft only as related work (Hilliard 1987; Ganie–Rather–Das 2023) were dropped as closed-access antecedent credits; nothing depends on them. |

## Scope of assurance

Internal AI-assisted review only: written proofs checked by AI reviewers, computations replayed, cited sources read.
No human specialist review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
Literature searches were of limited scope (MathSciNet and zbMATH were not searched; see
[RESEARCH_GATES.json](RESEARCH_GATES.json)).

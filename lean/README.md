# Invariance formalization

This Lean project formalizes the CLR invariance characterization from the
first section of the supplement.

The main tether declarations are:

```lean
Invariance.clr_uniqueness_theorem
Invariance.axiom_necessity_tether
Invariance.shifted_pflogpf_variance_scale
```

The first states that, in dimension `D + 2` (Lean's constructive form of
`D >= 2`), any rank-monotone additive map on log-coordinates that is
equivariant under coordinate relabeling, invariant on constant-depth shifts,
and calibrated on the first standard basis vector is exactly `Invariance.clr`.
The third proves the algebraic relationship between the software scale factor,
the count-scale pseudocount, and `alpha`. The second packages the five
dimension-3 independence witnesses for the LaTeX proposition that all five
axioms are necessary. The rank-monotonicity witness constructs a
non-order-preserving additive real map by extending `{1, sqrt 2}` to a
`ℚ`-basis of `ℝ`.

The proof includes the standard one-dimensional Cauchy regularity lemma that an
additive real function positive on the positive half-line is linear, recorded
as `Invariance.additive_positive_linear`.

Check the project with:

```bash
lake build
```

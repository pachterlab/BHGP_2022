# Invariance formalization

This Lean project formalizes the CLR invariance characterization from the
first section of the supplement.

The main result is:

```lean
Invariance.clr_characterization
```

It states that, in dimension `D + 2` (Lean's constructive form of
`D >= 2`), any rank-monotone additive map on log-coordinates that is
equivariant under coordinate relabeling, invariant on constant-depth shifts,
and calibrated on the first standard basis vector is exactly `Invariance.clr`.

The proof uses the standard one-dimensional Cauchy regularity lemma that an
additive real function positive on the positive half-line is linear, recorded
as `Invariance.additive_positive_linear`.

Check the project with:

```bash
lake build
```

import Mathlib

open scoped BigOperators

namespace Invariance

/-- Coordinate-wise permutation action on `Fin n -> R`. -/
def permute {n : Nat} (sigma : Equiv.Perm (Fin n)) (x : Fin n -> ℝ) : Fin n -> ℝ :=
  fun i => x (sigma i)

/-- The standard basis vector with a `1` in coordinate `i`. -/
def basis {n : Nat} (i : Fin n) : Fin n -> ℝ :=
  fun k => if k = i then 1 else 0

/-- The centered log-ratio projection on the log-domain. -/
noncomputable def clr {n : Nat} (x : Fin n -> ℝ) : Fin n -> ℝ :=
  x - ((∑ j : Fin n, x j) / (n : ℝ)) • (1 : Fin n -> ℝ)

/-- The first coordinate in a dimension known to be at least two. -/
def first (D : Nat) : Fin (D + 2) :=
  ⟨0, by omega⟩

/-- The second coordinate in a dimension known to be at least two. -/
def second (D : Nat) : Fin (D + 2) :=
  ⟨1, by omega⟩

lemma first_ne_second (D : Nat) : first D ≠ second D := by
  intro h
  have hv := congrArg Fin.val h
  simp [first, second] at hv

lemma sum_smul_basis {n : Nat} (x : Fin n -> ℝ) :
    (∑ j : Fin n, x j • basis j) = x := by
  classical
  ext k
  simp [basis]

lemma sum_basis_eq_one {n : Nat} :
    (∑ j : Fin n, basis j) = (1 : Fin n -> ℝ) := by
  classical
  ext k
  simp [basis]

lemma swap_apply_eq_left_iff {α : Type*} [DecidableEq α] (i j k : α) :
    Equiv.swap i j k = i ↔ k = j := by
  constructor
  · intro h
    have h' := congrArg (Equiv.swap i j) h
    simpa using h'
  · intro h
    subst h
    simp

lemma permute_swap_basis {n : Nat} (i j : Fin n) :
    permute (Equiv.swap i j) (basis i) = basis j := by
  classical
  ext k
  have hiff := swap_apply_eq_left_iff i j k
  by_cases hk : k = j
  · have hs : Equiv.swap i j k = i := hiff.mpr hk
    simp [permute, basis, hk]
  · have hs : Equiv.swap i j k ≠ i := by
      intro hs
      exact hk (hiff.mp hs)
    simp [permute, basis, hk, hs]

lemma permute_smul {n : Nat} (sigma : Equiv.Perm (Fin n)) (r : ℝ)
    (x : Fin n -> ℝ) :
    permute sigma (r • x) = r • permute sigma x := by
  ext k
  rfl

lemma permute_swap_basis_of_ne {n : Nat} {i p q : Fin n}
    (hp : i ≠ p) (hq : i ≠ q) :
    permute (Equiv.swap p q) (basis i) = basis i := by
  classical
  ext k
  have hfix : Equiv.swap p q i = i := Equiv.swap_apply_of_ne_of_ne hp hq
  by_cases hk : k = i
  · subst hk
    simp [permute, basis, hfix]
  · have hs : Equiv.swap p q k ≠ i := by
      intro hs
      have h' := congrArg (Equiv.swap p q) hs
      have hki : k = i := by
        simpa [hfix] using h'
      exact hk hki
    simp [permute, basis, hk, hs]

lemma sum_mul_if_eq_add {n : Nat} (x : Fin n -> ℝ) (i : Fin n) (a b : ℝ) :
    (∑ j : Fin n, x j * (if j = i then a else b)) =
      (a - b) * x i + b * (∑ j : Fin n, x j) := by
  classical
  let s : Finset (Fin n) := Finset.univ
  have hmul :
      Finset.sum (s := s.erase i) (fun j => x j * b) =
        b * Finset.sum (s := s.erase i) (fun j => x j) := by
    simpa [mul_comm, mul_left_comm, mul_assoc] using
      (Finset.sum_mul (s := s.erase i) (f := x) b).symm
  have hdecomp :
      (∑ j : Fin n, x j * (if j = i then a else b)) =
        x i * a + b * Finset.sum (s := s.erase i) (fun j => x j) := by
    calc
      (∑ j : Fin n, x j * (if j = i then a else b))
          = Finset.sum (s := s.erase i)
              (fun j => x j * (if j = i then a else b)) + x i * a := by
              simp [s, add_comm]
      _ = Finset.sum (s := s.erase i) (fun j => x j * b) + x i * a := by
            congr 1
            refine Finset.sum_congr rfl ?_
            intro j hj
            have hji : j ≠ i := (Finset.mem_erase.mp hj).1
            simp [hji]
      _ = b * Finset.sum (s := s.erase i) (fun j => x j) + x i * a := by
            rw [hmul]
      _ = x i * a + b * Finset.sum (s := s.erase i) (fun j => x j) := by
            ring
  have huniv :
      (∑ j : Fin n, x j) =
        x i + Finset.sum (s := s.erase i) (fun j => x j) := by
    have hi_mem : i ∈ s := by
      simp [s]
    have hsum :=
      Finset.sum_erase_add (s := s) (f := x) hi_mem
    linarith [hsum]
  calc
    (∑ j : Fin n, x j * (if j = i then a else b))
        = x i * a + b * Finset.sum (s := s.erase i) (fun j => x j) := hdecomp
    _ = (a - b) * x i + b *
          (x i + Finset.sum (s := s.erase i) (fun j => x j)) := by
          ring
    _ = (a - b) * x i + b * (∑ j : Fin n, x j) := by
          rw [huniv]

/--
Standard one-dimensional Cauchy regularity lemma used in the monotone CLR proof.
An additive real function that is positive on the positive half-line is linear.
-/
axiom additive_positive_linear (h : ℝ → ℝ)
    (hadd : ∀ x y : ℝ, h (x + y) = h x + h y)
    (hpos : ∀ x : ℝ, 0 < x → 0 < h x) :
    ∀ x : ℝ, h x = h 1 * x

/--
Additive characterization of CLR on the log-domain.

The formal dimension is `D + 2`, which is Lean's constructive way of saying
that the number of coordinates is at least two.  The result is the first
section of the supplement after passing to logarithmic coordinates: a
rank-monotone additive, permutation-equivariant, scale-invariant, naturally
calibrated map is the centered log-ratio projection.
-/
theorem clr_characterization (D : Nat)
    (T : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (hmono :
      ∀ x : Fin (D + 2) -> ℝ, ∀ i j : Fin (D + 2),
        x i > x j ↔ T x i > T x j)
    (hadd : ∀ x y : Fin (D + 2) -> ℝ, T (x + y) = T x + T y)
    (hperm :
      ∀ sigma : Equiv.Perm (Fin (D + 2)), ∀ x : Fin (D + 2) -> ℝ,
        T (permute sigma x) = permute sigma (T x))
    (hscale : ∀ c : ℝ, T (fun _ : Fin (D + 2) => c) = 0)
    (hcal : T (basis (first D)) (first D) =
      (D + 1 : ℝ) / (D + 2 : ℝ)) :
    ∀ x : Fin (D + 2) -> ℝ,
      T x = clr x := by
  classical
  let n := D + 2
  let i0 : Fin n := first D
  let i1 : Fin n := second D
  have hi01 : i0 ≠ i1 := by
    simpa [i0, i1] using first_ne_second D
  have hdim_ne : (D + 2 : ℝ) ≠ 0 := by positivity
  have hdim1_ne : (D + 1 : ℝ) ≠ 0 := by positivity

  let F : (Fin n -> ℝ) →+ (Fin n -> ℝ) := {
    toFun := T
    map_zero' := by
      ext i
      have hcoord : T (0 : Fin n -> ℝ) i =
          T (0 : Fin n -> ℝ) i + T (0 : Fin n -> ℝ) i := by
        simpa using congrArg (fun v => v i)
          (hadd (0 : Fin n -> ℝ) (0 : Fin n -> ℝ))
      have hzero : T (0 : Fin n -> ℝ) i = 0 := by
        linarith
      simpa using hzero
    map_add' := hadd
  }

  have hFT : ∀ x : Fin n -> ℝ, F x = T x := by
    intro x
    rfl
  have hpermF :
      ∀ sigma : Equiv.Perm (Fin n), ∀ x : Fin n -> ℝ,
        F (permute sigma x) = permute sigma (F x) := by
    intro sigma x
    simpa [F] using hperm sigma x

  let f1 : ℝ → ℝ := fun r => F (r • basis i0) i0
  let g1 : ℝ → ℝ := fun r => F (r • basis i0) i1
  let h1 : ℝ → ℝ := fun r => f1 r - g1 r

  have hfadd : ∀ x y : ℝ, f1 (x + y) = f1 x + f1 y := by
    intro x y
    have hs : (x + y) • basis i0 = x • basis i0 + y • basis i0 := by
      ext k
      by_cases hk : k = i0
      · simp [basis, hk]
      · simp [basis, hk]
    calc
      f1 (x + y) = F ((x + y) • basis i0) i0 := rfl
      _ = F (x • basis i0 + y • basis i0) i0 := by rw [hs]
      _ = (F (x • basis i0) + F (y • basis i0)) i0 := by
            rw [F.map_add]
      _ = f1 x + f1 y := rfl

  have hgadd : ∀ x y : ℝ, g1 (x + y) = g1 x + g1 y := by
    intro x y
    have hs : (x + y) • basis i0 = x • basis i0 + y • basis i0 := by
      ext k
      by_cases hk : k = i0
      · simp [basis, hk]
      · simp [basis, hk]
    calc
      g1 (x + y) = F ((x + y) • basis i0) i1 := rfl
      _ = F (x • basis i0 + y • basis i0) i1 := by rw [hs]
      _ = (F (x • basis i0) + F (y • basis i0)) i1 := by
            rw [F.map_add]
      _ = g1 x + g1 y := rfl

  have hhadd : ∀ x y : ℝ, h1 (x + y) = h1 x + h1 y := by
    intro x y
    dsimp [h1]
    rw [hfadd, hgadd]
    ring

  have hhpos : ∀ x : ℝ, 0 < x → 0 < h1 x := by
    intro x hx
    have hb0 : (x • basis i0) i0 > (x • basis i0) i1 := by
      simp [basis, hi01.symm, hx]
    have hm := (hmono (x • basis i0) i0 i1).mp hb0
    simpa [h1, f1, g1, F] using sub_pos.mpr hm

  have hhlin : ∀ x : ℝ, h1 x = h1 1 * x :=
    additive_positive_linear h1 hhadd hhpos

  let alpha : ℝ := h1 1

  have halpha_pos : 0 < alpha := by
    simpa [alpha] using hhpos 1 zero_lt_one

  have hEi : ∀ i : Fin n,
      F (basis i) = permute (Equiv.swap i0 i) (F (basis i0)) := by
    intro i
    have h := hpermF (Equiv.swap i0 i) (basis i0)
    simpa [permute_swap_basis] using h

  have hoff0 : ∀ r : ℝ, ∀ q : Fin n, q ≠ i0 ->
      F (r • basis i0) q = g1 r := by
    intro r q hq
    by_cases hqi : q = i1
    · subst hqi
      rfl
    · have hfix :
          permute (Equiv.swap i1 q) (r • basis i0) = r • basis i0 := by
        rw [permute_smul, permute_swap_basis_of_ne hi01 hq.symm]
      have hp := hpermF (Equiv.swap i1 q) (r • basis i0)
      have hp0 : F (r • basis i0) =
          permute (Equiv.swap i1 q) (F (r • basis i0)) := by
        simpa [hfix] using hp
      have hc := congrArg (fun v => v i1) hp0
      have hb' : g1 r = F (r • basis i0) q := by
        simpa [g1, permute, hqi] using hc
      exact hb'.symm

  have hscale_fg : ∀ r : ℝ, f1 r + (D + 1 : ℝ) * g1 r = 0 := by
    intro r
    have hone : (r • (1 : Fin n -> ℝ)) = (fun _ : Fin n => r) := by
      ext k
      simp
    have hsum_basis : (r • (1 : Fin n -> ℝ)) = ∑ j : Fin n, r • basis j := by
      rw [← sum_basis_eq_one (n := n)]
      simp [Finset.smul_sum]
    have hzero : F (r • (1 : Fin n -> ℝ)) = 0 := by
      simpa [F, hone] using hscale r
    have hsum :
        F (r • (1 : Fin n -> ℝ)) i0 =
          (∑ j : Fin n, F (r • basis j) i0) := by
      rw [hsum_basis]
      simpa using congrArg (fun v => v i0)
        (map_sum F (fun j : Fin n => r • basis j) (Finset.univ : Finset (Fin n)))
    have hcoeff :
        (∑ j : Fin n, F (r • basis j) i0) =
          (∑ j : Fin n, (if j = i0 then f1 r else g1 r)) := by
      refine Finset.sum_congr rfl ?_
      intro j _
      by_cases hji : j = i0
      · subst hji
        simp [f1]
      · have h := hpermF (Equiv.swap i0 j) (r • basis i0)
        have hbasis : permute (Equiv.swap i0 j) (r • basis i0) = r • basis j := by
          rw [permute_smul, permute_swap_basis]
        have hcoord := congrArg (fun v => v i0) h
        have hswap : Equiv.swap i0 j i0 = j := by simp
        have hval : F (r • basis j) i0 = g1 r := by
          have htmp : F (r • basis j) i0 = F (r • basis i0) j := by
            simpa [hbasis, permute, hswap] using hcoord
          rw [htmp]
          exact hoff0 r j hji
        simp [hji, hval]
    have hones : (∑ j : Fin n, (1 : ℝ)) = (D + 2 : ℝ) := by
      simp [n]
    have hsum_if :
        (∑ j : Fin n, (if j = i0 then f1 r else g1 r)) =
          f1 r + (D + 1 : ℝ) * g1 r := by
      calc
        (∑ j : Fin n, (if j = i0 then f1 r else g1 r))
            = (∑ j : Fin n, (1 : ℝ) * (if j = i0 then f1 r else g1 r)) := by
                simp
        _ = (f1 r - g1 r) * (1 : ℝ) + g1 r * (∑ j : Fin n, (1 : ℝ)) := by
              simpa using
                sum_mul_if_eq_add (fun _ : Fin n => (1 : ℝ)) i0 (f1 r) (g1 r)
        _ = f1 r + (D + 1 : ℝ) * g1 r := by
              rw [hones]
              ring
    have hleft : F (r • (1 : Fin n -> ℝ)) i0 = 0 := by
      simpa using congrArg (fun v => v i0) hzero
    linarith [hsum, hcoeff, hsum_if, hleft]

  have hf_formula : ∀ r : ℝ, f1 r = alpha * ((D + 1 : ℝ) / (D + 2 : ℝ)) * r := by
    intro r
    have hdiff : f1 r - g1 r = alpha * r := by
      simpa [h1, alpha] using hhlin r
    have hsum := hscale_fg r
    have hD : (D + 2 : ℝ) ≠ 0 := by positivity
    field_simp [hD]
    nlinarith [hdiff, hsum]

  have hg_formula : ∀ r : ℝ, g1 r = -alpha * (1 / (D + 2 : ℝ)) * r := by
    intro r
    have hdiff : f1 r - g1 r = alpha * r := by
      simpa [h1, alpha] using hhlin r
    have hsum := hscale_fg r
    have hD : (D + 2 : ℝ) ≠ 0 := by positivity
    field_simp [hD]
    nlinarith [hdiff, hsum]

  have hcalF : f1 1 = (D + 1 : ℝ) / (D + 2 : ℝ) := by
    simpa [f1, F, i0, n] using hcal

  have halpha : alpha = 1 := by
    have hf1 := hf_formula 1
    have hf1mul : f1 1 * (D + 2 : ℝ) = alpha * (D + 1 : ℝ) := by
      field_simp [hdim_ne] at hf1
      nlinarith [hf1]
    have hcalFmul : f1 1 * (D + 2 : ℝ) = (D + 1 : ℝ) := by
      field_simp [hdim_ne] at hcalF
      nlinarith [hcalF]
    have hdim1_pos : 0 < (D + 1 : ℝ) := by positivity
    nlinarith [hf1mul, hcalFmul, hdim1_pos]

  have hcoord_basis : ∀ r : ℝ, ∀ j i : Fin n,
      F (r • basis j) i =
        r * (if j = i then (D + 1 : ℝ) / (D + 2 : ℝ)
             else -(1 / (D + 2 : ℝ))) := by
    intro r j i
    by_cases hji : j = i
    · subst j
      have h := hpermF (Equiv.swap i0 i) (r • basis i0)
      have hbasis : permute (Equiv.swap i0 i) (r • basis i0) = r • basis i := by
        rw [permute_smul, permute_swap_basis]
      have hcoord := congrArg (fun v => v i) h
      have hswap : Equiv.swap i0 i i = i0 := by
        exact (swap_apply_eq_left_iff i0 i i).mpr rfl
      calc
        F (r • basis i) i = F (r • basis i0) i0 := by
          simpa [hbasis, permute, hswap] using hcoord
        _ = r * ((D + 1 : ℝ) / (D + 2 : ℝ)) := by
          change f1 r = r * ((D + 1 : ℝ) / (D + 2 : ℝ))
          rw [hf_formula, halpha]
          ring
        _ = r * (if i = i then (D + 1 : ℝ) / (D + 2 : ℝ)
             else -(1 / (D + 2 : ℝ))) := by simp
    · have h := hpermF (Equiv.swap i0 j) (r • basis i0)
      have hbasis : permute (Equiv.swap i0 j) (r • basis i0) = r • basis j := by
        rw [permute_smul, permute_swap_basis]
      have hcoord := congrArg (fun v => v i) h
      have hswap_ne : Equiv.swap i0 j i ≠ i0 := by
        intro hs
        have h' := congrArg (Equiv.swap i0 j) hs
        have hswap0 : Equiv.swap i0 j i0 = j := by simp
        have hij : i = j := by
          simpa [hswap0] using h'
        exact hji hij.symm
      have hoff0' : F (r • basis i0) (Equiv.swap i0 j i) = g1 r :=
        hoff0 r (Equiv.swap i0 j i) hswap_ne
      calc
        F (r • basis j) i =
            F (r • basis i0) (Equiv.swap i0 j i) := by
              simpa [hbasis, permute] using hcoord
        _ = r * (-(1 / (D + 2 : ℝ))) := by
              rw [hoff0', hg_formula, halpha]
              ring
        _ = r * (if j = i then (D + 1 : ℝ) / (D + 2 : ℝ)
             else -(1 / (D + 2 : ℝ))) := by simp [hji]

  intro x
  ext i
  calc
    T x i = F x i := by rfl
    _ = F (∑ j : Fin n, x j • basis j) i := by
          rw [sum_smul_basis]
    _ = (∑ j : Fin n, F (x j • basis j) i) := by
          simpa using congrArg (fun v => v i)
            (map_sum F (fun j : Fin n => x j • basis j)
              (Finset.univ : Finset (Fin n)))
    _ = (∑ j : Fin n,
          x j * (if j = i then (D + 1 : ℝ) / (D + 2 : ℝ)
             else -(1 / (D + 2 : ℝ)))) := by
          refine Finset.sum_congr rfl ?_
          intro j _
          exact hcoord_basis (x j) j i
    _ = (((D + 1 : ℝ) / (D + 2 : ℝ)) - (-(1 / (D + 2 : ℝ)))) * x i
          + (-(1 / (D + 2 : ℝ))) * (∑ j : Fin n, x j) := by
          simpa using
            sum_mul_if_eq_add x i ((D + 1 : ℝ) / (D + 2 : ℝ)) (-(1 / (D + 2 : ℝ)))
    _ = x i - (∑ j : Fin n, x j) / (D + 2 : ℝ) := by
          field_simp [hdim_ne]
          ring
    _ = clr x i := by
          simp [clr, n]

end Invariance

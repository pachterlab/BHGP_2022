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

/-- Graph label: the first and second coordinates are distinct in dimension `D + 2`. -/
lemma first_ne_second (D : Nat) : first D ≠ second D := by
  intro h
  have hv := congrArg Fin.val h
  simp [first, second] at hv

/-- Graph label: finite-coordinate expansion in the standard basis. -/
lemma sum_smul_basis {n : Nat} (x : Fin n -> ℝ) :
    (∑ j : Fin n, x j • basis j) = x := by
  classical
  ext k
  simp [basis]

/-- Graph label: the sum of standard basis vectors is the all-ones vector. -/
lemma sum_basis_eq_one {n : Nat} :
    (∑ j : Fin n, basis j) = (1 : Fin n -> ℝ) := by
  classical
  ext k
  simp [basis]

/-- Graph label: coordinate criterion for a transposition hitting its left endpoint. -/
lemma swap_apply_eq_left_iff {α : Type*} [DecidableEq α] (i j k : α) :
    Equiv.swap i j k = i ↔ k = j := by
  constructor
  · intro h
    have h' := congrArg (Equiv.swap i j) h
    simpa using h'
  · intro h
    subst h
    simp

/-- Graph label: swapping coordinates sends one standard basis vector to the other. -/
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

/-- Graph label: coordinate permutations commute with scalar multiplication. -/
lemma permute_smul {n : Nat} (sigma : Equiv.Perm (Fin n)) (r : ℝ)
    (x : Fin n -> ℝ) :
    permute sigma (r • x) = r • permute sigma x := by
  ext k
  rfl

/-- Graph label: swapping two other coordinates fixes a standard basis vector. -/
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

/-- Graph label: closed form for a finite sum with one distinguished coordinate. -/
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
Graph label: Cauchy regularity for additive monotone real maps.

Standard one-dimensional Cauchy regularity lemma used in the monotone CLR proof.
An additive real function that is positive on the positive half-line is linear.
-/
theorem additive_positive_linear (h : ℝ → ℝ)
    (hadd : ∀ x y : ℝ, h (x + y) = h x + h y)
    (hpos : ∀ x : ℝ, 0 < x → 0 < h x) :
    ∀ x : ℝ, h x = h 1 * x := by
  let H : ℝ →+ ℝ := {
    toFun := h
    map_zero' := by
      have hz : h 0 = h 0 + h 0 := by simpa using hadd 0 0
      linarith
    map_add' := hadd
  }
  have h1pos : 0 < h 1 := hpos 1 zero_lt_one
  have hmono : Monotone h := by
    intro x y hxy
    rcases lt_or_eq_of_le hxy with hlt | rfl
    · have hd : 0 < y - x := sub_pos.mpr hlt
      have hp : 0 < h (y - x) := hpos (y - x) hd
      have hy : h y = h x + h (y - x) := by
        calc
          h y = h (x + (y - x)) := by ring_nf
          _ = h x + h (y - x) := hadd x (y - x)
      linarith
    · exact le_rfl
  have hrat : ∀ q : ℚ, h (q : ℝ) = h 1 * (q : ℝ) := by
    intro q
    have hm := map_rat_smul H q (1 : ℝ)
    simpa [H, Rat.smul_def, mul_comm] using hm
  intro x
  apply le_antisymm
  · by_contra hnot
    have hlt : h 1 * x < h x := lt_of_not_ge hnot
    have hxlt : x < h x / h 1 := by
      rw [lt_div_iff₀ h1pos]
      simpa [mul_comm] using hlt
    obtain ⟨q, hxq, hq⟩ := exists_rat_btwn hxlt
    have hle : h x ≤ h (q : ℝ) := hmono hxq.le
    have hq_lt : h 1 * (q : ℝ) < h x := by
      have hmul := mul_lt_mul_of_pos_left hq h1pos
      field_simp [ne_of_gt h1pos] at hmul
      simpa [mul_comm, mul_left_comm, mul_assoc] using hmul
    rw [hrat q] at hle
    linarith
  · by_contra hnot
    have hlt : h x < h 1 * x := lt_of_not_ge hnot
    have hdivlt : h x / h 1 < x := by
      rw [div_lt_iff₀ h1pos]
      simpa [mul_comm] using hlt
    obtain ⟨q, hq, hqx⟩ := exists_rat_btwn hdivlt
    have hle : h (q : ℝ) ≤ h x := hmono hqx.le
    have hx_lt_q : h x < h 1 * (q : ℝ) := by
      have hmul := mul_lt_mul_of_pos_left hq h1pos
      field_simp [ne_of_gt h1pos] at hmul
      simpa [mul_comm, mul_left_comm, mul_assoc] using hmul
    rw [hrat q] at hle
    linarith

/--
Graph label: main CLR characterization.

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

/--
Graph label: LaTeX tether for CLR uniqueness.

Named tether for Supplementary Note Theorem `thm:clr-uniqueness`.

This is the same formal statement as `clr_characterization`: after passing to
log-coordinates, the five axioms force the centered log-ratio projection.
-/
theorem clr_uniqueness_theorem (D : Nat)
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
      T x = clr x :=
  clr_characterization D T hmono hadd hperm hscale hcal

/-- The log1pPF expression with software scale factor `K`. -/
noncomputable def log1pPFScale (K s : ℝ) (x : ℝ) : ℝ :=
  Real.log (1 + K * (x / s))

/-- The same expression written with the corresponding count-scale pseudocount. -/
noncomputable def shiftedLogCountScale (K s : ℝ) (x : ℝ) : ℝ :=
  Real.log ((K / s) * (x + s / K))

/--
Graph label: scale-factor/count-pseudocount identity.

The algebra behind the count-scale pseudocount: applying a software scale
factor `K` after depth normalization is the same, up to the multiplicative
constant `K/s` inside the logarithm, as adding the count-scale pseudocount
`s/K` before taking the logarithm.
-/
theorem scale_factor_count_pseudocount_identity
    {K s x : ℝ} (hK : K ≠ 0) (hs : s ≠ 0) :
    1 + K * (x / s) = (K / s) * (x + s / K) := by
  field_simp [hK, hs]
  ring

/-- Graph label: logarithmic form of `scale_factor_count_pseudocount_identity`. -/
theorem log1pPFScale_eq_shiftedLogCountScale
    {K s x : ℝ} (hK : K ≠ 0) (hs : s ≠ 0) :
    log1pPFScale K s x = shiftedLogCountScale K s x := by
  unfold log1pPFScale shiftedLogCountScale
  rw [scale_factor_count_pseudocount_identity (K := K) (s := s) (x := x) hK hs]

/--
Graph label: solve the software scale factor from `alpha`.

If the count-scale pseudocount is `1 / (4 * alpha)`, then the software scale
factor is `K = 4 * alpha * s`.
-/
theorem scale_factor_from_alpha
    {K s alpha : ℝ} (hK : K ≠ 0) (halpha : alpha ≠ 0)
    (h : s / K = 1 / (4 * alpha)) :
    K = 4 * alpha * s := by
  have h4a : 4 * alpha ≠ 0 := by
    exact mul_ne_zero (by norm_num) halpha
  field_simp [hK, h4a] at h
  nlinarith

/--
Graph label: LaTeX tether for shifted PFlog scale.

Named tether for Supplementary Note Proposition `prop:shifted-pflog-scale`.
It packages the two exact algebraic facts used in the proposition.
-/
theorem shifted_pflog_variance_scale
    {K s alpha x : ℝ} (hK : K ≠ 0) (hs : s ≠ 0) (halpha : alpha ≠ 0)
    (h : s / K = 1 / (4 * alpha)) :
    log1pPFScale K s x = shiftedLogCountScale K s x ∧
      K = 4 * alpha * s := by
  constructor
  · exact log1pPFScale_eq_shiftedLogCountScale (K := K) (s := s) (x := x) hK hs
  · exact scale_factor_from_alpha (K := K) (s := s) (alpha := alpha) hK halpha h

/-- Coordinatewise logarithm of a positive-domain vector. -/
noncomputable def logCoords {n : Nat} (z : Fin n -> ℝ) : Fin n -> ℝ :=
  fun i => Real.log (z i)

/-- CLR applied after coordinatewise logarithm on the positive domain. -/
noncomputable def positiveCLR {n : Nat} (z : Fin n -> ℝ) : Fin n -> ℝ :=
  clr (logCoords z)

/-- Add a fixed normalized shift to every coordinate. -/
def shiftedComposition {n : Nat} (tau : ℝ) (u : Fin n -> ℝ) : Fin n -> ℝ :=
  fun i => u i + tau

/-- The shifted CLR transform with fixed normalized shift `tau`. -/
noncomputable def shiftedCLR {n : Nat} (tau : ℝ) (u : Fin n -> ℝ) : Fin n -> ℝ :=
  positiveCLR (shiftedComposition tau u)

/-- The positive-domain condition for the shifted composition. -/
def ShiftedPositive {n : Nat} (tau : ℝ) (u : Fin n -> ℝ) : Prop :=
  ∀ i : Fin n, 0 < u i + tau

/--
Graph label: shifted-domain pullback of CLR uniqueness.

If a positive-domain transform `S` is obtained by pulling back a log-domain
transform `T` through coordinatewise logarithms, and `T` satisfies the CLR
axioms in log-coordinates, then `S` is the ordinary CLR of the logged positive
vector. This is the formal restriction step used in Supplementary Note
Proposition `prop:shifted-pflogpf-scale`.
-/
theorem positive_domain_pullback_clr (D : Nat)
    (S : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (T : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (hS : ∀ z : Fin (D + 2) -> ℝ,
      (∀ i : Fin (D + 2), 0 < z i) -> S z = T (logCoords z))
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
    ∀ z : Fin (D + 2) -> ℝ, (∀ i : Fin (D + 2), 0 < z i) ->
      S z = positiveCLR z := by
  intro z hz
  rw [hS z hz]
  unfold positiveCLR
  exact clr_uniqueness_theorem D T hmono hadd hperm hscale hcal (logCoords z)

/--
Graph label: shifted-composition pullback of CLR uniqueness.

For a fixed normalized shift `tau`, if a transform on unshifted compositions is
the restriction of a positive-domain transform satisfying the preceding
pullback statement after applying `u ↦ u + tau`, then it is exactly shifted
CLR. This is the Lean formalization of the shifted-domain pullback statement.
-/
theorem shifted_domain_pullback_clr (D : Nat) (tau : ℝ)
    (S : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (T : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (hS : ∀ u : Fin (D + 2) -> ℝ, ShiftedPositive tau u ->
      S u = T (logCoords (shiftedComposition tau u)))
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
    ∀ u : Fin (D + 2) -> ℝ, ShiftedPositive tau u ->
      S u = shiftedCLR tau u := by
  intro u hu
  rw [hS u hu]
  unfold shiftedCLR positiveCLR
  exact clr_uniqueness_theorem D T hmono hadd hperm hscale hcal
    (logCoords (shiftedComposition tau u))

/-- Rank monotonicity on the log-coordinate domain. -/
def RankMonotone {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∀ x : Fin n -> ℝ, ∀ i j : Fin n,
    x i > x j ↔ T x i > T x j

/-- Perturbation additivity after passing to log-coordinates. -/
def PerturbationAdditive {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∀ x y : Fin n -> ℝ, T (x + y) = T x + T y

/-- Equivariance under relabeling coordinates. -/
def RelabelingEquivariant {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∀ sigma : Equiv.Perm (Fin n), ∀ x : Fin n -> ℝ,
    T (permute sigma x) = permute sigma (T x)

/-- Invariance under adding a constant depth shift in log-coordinates. -/
def ScaleInvariant {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∀ c : ℝ, T (fun _ : Fin n => c) = 0

/-- Full scale invariance in log-coordinates, not using additivity. -/
def ScaleInvariantFull {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∀ x : Fin n -> ℝ, ∀ c : ℝ, T (x + fun _ : Fin n => c) = T x

/-- Natural-log calibration in dimension `D + 2`. -/
def NaturallyCalibrated (D : Nat)
    (T : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ)) : Prop :=
  T (basis (first D)) (first D) = (D + 1 : ℝ) / (D + 2 : ℝ)

/-- The transform is not the CLR projection. -/
def DiffersFromCLR {n : Nat}
    (T : (Fin n -> ℝ) -> (Fin n -> ℝ)) : Prop :=
  ∃ x : Fin n -> ℝ, T x ≠ clr x

/-- A positive scalar multiple of CLR is rank-monotone in dimension 3. -/
lemma positive_smul_clr_rank_monotone_fin3 {a : ℝ} (ha : 0 < a) :
    RankMonotone (n := 3) (fun x => a • clr x) := by
  intro x i j
  have hiff : clr x i > clr x j ↔ x i > x j := by
    constructor
    · intro h
      simp [clr] at h
      linarith
    · intro h
      simp [clr]
      linarith
  constructor
  · intro h
    have hc : clr x i > clr x j := hiff.mpr h
    simpa [Pi.smul_apply] using mul_lt_mul_of_pos_left hc ha
  · intro h
    have hc : clr x i > clr x j := by
      have h' : a * clr x i > a * clr x j := by
        simpa [Pi.smul_apply] using h
      nlinarith
    exact hiff.mp hc

/-- A scalar multiple of CLR is additive in dimension 3. -/
lemma smul_clr_additive_fin3 (a : ℝ) :
    PerturbationAdditive (n := 3) (fun x => a • clr x) := by
  intro x y
  ext i
  fin_cases i
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]
    ring
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]
    ring
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]
    ring

/-- A scalar multiple of CLR is equivariant in dimension 3. -/
lemma smul_clr_equivariant_fin3 (a : ℝ) :
    RelabelingEquivariant (n := 3) (fun x => a • clr x) := by
  intro sigma x
  ext i
  have hsum :
      (∑ j : Fin 3, x (sigma j)) = ∑ j : Fin 3, x j := by
    exact Equiv.sum_comp sigma x
  simp [clr, permute, Pi.smul_apply, hsum]

/-- A scalar multiple of CLR is scale-invariant in dimension 3. -/
lemma smul_clr_scale_invariant_fin3 (a : ℝ) :
    ScaleInvariant (n := 3) (fun x => a • clr x) := by
  intro c
  ext i
  fin_cases i
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]
  · norm_num [clr, Pi.smul_apply, Fin.sum_univ_three]

/-- Graph label: independence witness when calibration is omitted. -/
theorem calibration_axiom_necessary :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ RelabelingEquivariant T ∧
        ScaleInvariant T ∧ DiffersFromCLR T := by
  let T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ) := fun x => (2 : ℝ) • clr x
  refine ⟨T, ?_, ?_, ?_, ?_, ?_⟩
  · exact positive_smul_clr_rank_monotone_fin3 (by norm_num)
  · exact smul_clr_additive_fin3 2
  · exact smul_clr_equivariant_fin3 2
  · exact smul_clr_scale_invariant_fin3 2
  · refine ⟨basis (⟨0, by norm_num⟩ : Fin 3), ?_⟩
    intro h
    have hc := congrArg (fun v : Fin 3 -> ℝ => v (⟨0, by norm_num⟩ : Fin 3)) h
    norm_num [T, clr, basis, Fin.sum_univ_three] at hc
    have h20 : (2 : Fin 3) ≠ 0 := by
      decide
    simp [h20] at hc
    norm_num at hc

/-- A positive scalar multiple of the identity is rank-monotone. -/
lemma positive_smul_id_rank_monotone_fin3 {a : ℝ} (ha : 0 < a) :
    RankMonotone (n := 3) (fun x : Fin 3 -> ℝ => a • x) := by
  intro x i j
  constructor
  · intro h
    simpa [Pi.smul_apply] using mul_lt_mul_of_pos_left h ha
  · intro h
    have h' : a * x i > a * x j := by
      simpa [Pi.smul_apply] using h
    nlinarith

/-- A scalar multiple of the identity is additive. -/
lemma smul_id_additive_fin3 (a : ℝ) :
    PerturbationAdditive (n := 3) (fun x : Fin 3 -> ℝ => a • x) := by
  intro x y
  ext i
  simp [Pi.smul_apply]
  ring

/-- A scalar multiple of the identity is equivariant. -/
lemma smul_id_equivariant_fin3 (a : ℝ) :
    RelabelingEquivariant (n := 3) (fun x : Fin 3 -> ℝ => a • x) := by
  intro sigma x
  ext i
  rfl

/-- Graph label: independence witness when scale invariance is omitted. -/
theorem scale_invariance_axiom_necessary :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ RelabelingEquivariant T ∧
        NaturallyCalibrated 1 T ∧ DiffersFromCLR T := by
  let T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ) := fun x => ((2 : ℝ) / 3) • x
  refine ⟨T, ?_, ?_, ?_, ?_, ?_⟩
  · exact positive_smul_id_rank_monotone_fin3 (by norm_num)
  · exact smul_id_additive_fin3 ((2 : ℝ) / 3)
  · exact smul_id_equivariant_fin3 ((2 : ℝ) / 3)
  · norm_num [NaturallyCalibrated, T, basis, first]
  · refine ⟨(fun _ : Fin 3 => (1 : ℝ)), ?_⟩
    intro h
    have hc := congrArg (fun v : Fin 3 -> ℝ => v (⟨0, by norm_num⟩ : Fin 3)) h
    norm_num [T, clr, Fin.sum_univ_three] at hc

/-- A non-uniform weighted centering used when relabeling equivariance is omitted. -/
noncomputable def weightedCenter3 (x : Fin 3 -> ℝ) : Fin 3 -> ℝ :=
  fun i => x i -
    ((1 : ℝ) / 3 * x (⟨0, by norm_num⟩ : Fin 3) +
      (1 : ℝ) / 2 * x (⟨1, by norm_num⟩ : Fin 3) +
      (1 : ℝ) / 6 * x (⟨2, by norm_num⟩ : Fin 3))

/-- Graph label: nonuniform weighted centering preserves coordinate order. -/
lemma weightedCenter3_rank_monotone :
    RankMonotone weightedCenter3 := by
  intro x i j
  constructor
  · intro h
    simp [weightedCenter3]
    linarith
  · intro h
    simp [weightedCenter3] at h
    linarith

/-- Graph label: nonuniform weighted centering is additive. -/
lemma weightedCenter3_additive :
    PerturbationAdditive weightedCenter3 := by
  intro x y
  ext i
  fin_cases i
  · norm_num [weightedCenter3]
    ring
  · norm_num [weightedCenter3]
    ring
  · norm_num [weightedCenter3]
    ring

/-- Graph label: nonuniform weighted centering kills constant depth shifts. -/
lemma weightedCenter3_scale_invariant :
    ScaleInvariant weightedCenter3 := by
  intro c
  ext i
  fin_cases i
  · norm_num [weightedCenter3]
    ring
  · norm_num [weightedCenter3]
    ring
  · norm_num [weightedCenter3]
    ring

/-- Graph label: nonuniform weighted centering has the natural calibration. -/
lemma weightedCenter3_calibrated :
    NaturallyCalibrated 1 weightedCenter3 := by
  norm_num [NaturallyCalibrated, weightedCenter3, basis, first]

/-- Graph label: independence witness when relabeling equivariance is omitted. -/
theorem relabeling_equivariance_axiom_necessary :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ ScaleInvariant T ∧
        NaturallyCalibrated 1 T ∧ DiffersFromCLR T := by
  refine ⟨weightedCenter3, weightedCenter3_rank_monotone,
    weightedCenter3_additive, weightedCenter3_scale_invariant,
    weightedCenter3_calibrated, ?_⟩
  refine ⟨basis (⟨1, by norm_num⟩ : Fin 3), ?_⟩
  intro h
  have hc := congrArg (fun v : Fin 3 -> ℝ => v (⟨0, by norm_num⟩ : Fin 3)) h
  have h10 : (1 : Fin 3) ≠ 0 := by decide
  have h20 : (2 : Fin 3) ≠ 0 := by decide
  have h21 : (2 : Fin 3) ≠ 1 := by decide
  norm_num [weightedCenter3, clr, basis, Fin.sum_univ_three, h10, h20, h21] at hc

/-- Centering on the log-coordinate domain in dimension 3. -/
noncomputable def center3 (x : Fin 3 -> ℝ) : Fin 3 -> ℝ :=
  fun i => x i - (∑ j : Fin 3, x j) / 3

/-- A nonadditive, rank-monotone, scale-invariant calibrated transform. -/
noncomputable def cubicCentered3 (x : Fin 3 -> ℝ) : Fin 3 -> ℝ :=
  fun i => ((9 : ℝ) / 13) * (center3 x i + center3 x i ^ 3)

/-- Graph label: the centered cubic counterexample preserves coordinate order. -/
lemma cubicCentered3_rank_monotone :
    RankMonotone cubicCentered3 := by
  intro x i j
  let A : ℝ := (9 : ℝ) / 13
  have hA : 0 < A := by
    dsimp [A]
    positivity
  have hq : StrictMono (fun t : ℝ => t + t ^ 3) := by
    intro a b hab
    have hpow : a ^ 3 ≤ b ^ 3 :=
      (show Odd (3 : Nat) by norm_num).strictMono_pow.monotone hab.le
    linarith
  have hciff : center3 x i > center3 x j ↔ x i > x j := by
    constructor <;> intro h <;> simp [center3] at h ⊢ <;> linarith
  constructor
  · intro h
    have hc : center3 x i > center3 x j := hciff.mpr h
    have hq' :
        center3 x i + center3 x i ^ 3 >
          center3 x j + center3 x j ^ 3 :=
      hq hc
    simpa [cubicCentered3, A] using mul_lt_mul_of_pos_left hq' hA
  · intro h
    have hq' :
        center3 x i + center3 x i ^ 3 >
          center3 x j + center3 x j ^ 3 := by
      by_contra hnot
      have hle :
          center3 x i + center3 x i ^ 3 ≤
            center3 x j + center3 x j ^ 3 := le_of_not_gt hnot
      have hAle :
          A * (center3 x i + center3 x i ^ 3) ≤
            A * (center3 x j + center3 x j ^ 3) :=
        mul_le_mul_of_nonneg_left hle hA.le
      have hlt :
          A * (center3 x i + center3 x i ^ 3) >
            A * (center3 x j + center3 x j ^ 3) := by
        simpa [cubicCentered3, A] using h
      exact not_lt_of_ge hAle hlt
    exact hciff.mp (hq.lt_iff_lt.mp hq')

/-- Graph label: the centered cubic counterexample is relabeling equivariant. -/
lemma cubicCentered3_equivariant :
    RelabelingEquivariant cubicCentered3 := by
  intro sigma x
  ext i
  have hsum :
      (∑ j : Fin 3, x (sigma j)) = ∑ j : Fin 3, x j := by
    exact Equiv.sum_comp sigma x
  simp [cubicCentered3, center3, permute, hsum]

/-- Graph label: the centered cubic counterexample is fully scale invariant. -/
lemma cubicCentered3_scale_invariant_full :
    ScaleInvariantFull cubicCentered3 := by
  intro x c
  ext i
  have hsum :
      (∑ j : Fin 3, (x j + c)) = (∑ j : Fin 3, x j) + 3 * c := by
    norm_num [Fin.sum_univ_three]
    ring
  simp [cubicCentered3, center3, hsum]
  ring_nf

/-- Graph label: the centered cubic counterexample has the natural calibration. -/
lemma cubicCentered3_calibrated :
    NaturallyCalibrated 1 cubicCentered3 := by
  have h20 : (2 : Fin 3) ≠ 0 := by decide
  norm_num [NaturallyCalibrated, cubicCentered3, center3, basis, first,
    Fin.sum_univ_three, h20]

/-- Graph label: the centered cubic counterexample violates perturbation additivity. -/
lemma cubicCentered3_not_additive :
    ¬ PerturbationAdditive cubicCentered3 := by
  intro hadd
  have h := congrArg
    (fun v : Fin 3 -> ℝ => v (⟨0, by norm_num⟩ : Fin 3))
    (hadd (basis (⟨0, by norm_num⟩ : Fin 3))
      (basis (⟨0, by norm_num⟩ : Fin 3)))
  have h10 : (1 : Fin 3) ≠ 0 := by decide
  have h20 : (2 : Fin 3) ≠ 0 := by decide
  norm_num [cubicCentered3, center3, basis, Fin.sum_univ_three, h10, h20] at h

/-- Graph label: independence witness when perturbation additivity is omitted. -/
theorem perturbation_additivity_axiom_necessary :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ RelabelingEquivariant T ∧ ScaleInvariantFull T ∧
        NaturallyCalibrated 1 T ∧ ¬ PerturbationAdditive T := by
  exact ⟨cubicCentered3, cubicCentered3_rank_monotone, cubicCentered3_equivariant,
    cubicCentered3_scale_invariant_full, cubicCentered3_calibrated,
    cubicCentered3_not_additive⟩

/-- No rational multiple of `1` is `sqrt 2`. -/
lemma sqrt_two_not_rat_smul_one (q : ℚ) :
    q • (1 : ℝ) ≠ Real.sqrt 2 := by
  intro h
  have hq : (q : ℝ) = Real.sqrt 2 := by simpa using h
  exact irrational_sqrt_two.ne_rat q hq.symm

/-- The vectors `1` and `sqrt 2` are linearly independent over `ℚ`. -/
lemma sqrt_two_pair_linearIndependent :
    LinearIndependent ℚ ![(1 : ℝ), Real.sqrt 2] := by
  rw [LinearIndependent.pair_iff' (x := (1 : ℝ)) (y := Real.sqrt 2)]
  · intro q
    exact sqrt_two_not_rat_smul_one q
  · norm_num

/--
Graph label: construct the non-order-preserving additive real map.

A Cauchy-pathology witness used for the independence of rank monotonicity.
The construction extends `{1, sqrt 2}` to a `ℚ`-basis of `ℝ`, fixes all basis
directions except `sqrt 2`, and sends `sqrt 2` to `-sqrt 2`.
-/
theorem exists_rank_bad_additive_calibrated_phi :
    ∃ phi : ℝ -> ℝ,
      (∀ x y : ℝ, phi (x + y) = phi x + phi y) ∧
        phi ((2 : ℝ) / 3) = (2 : ℝ) / 3 ∧
          ∃ a b : ℝ, a > b ∧ ¬ phi a > phi b := by
  let v : Fin 2 -> ℝ := ![(1 : ℝ), Real.sqrt 2]
  have hv : LinearIndependent ℚ v := by
    simpa [v] using sqrt_two_pair_linearIndependent
  let s : Set ℝ := Set.range v
  have hs : LinearIndepOn ℚ id s := hv.linearIndepOn_id
  let B := Module.Basis.extend hs
  have h1mem : (1 : ℝ) ∈ s := by
    refine ⟨0, ?_⟩
    simp [v]
  have h2mem : Real.sqrt 2 ∈ s := by
    refine ⟨1, ?_⟩
    simp [v]
  let i1 : hs.extend (Set.subset_univ s) :=
    ⟨1, hs.subset_extend (Set.subset_univ s) h1mem⟩
  let i2 : hs.extend (Set.subset_univ s) :=
    ⟨Real.sqrt 2, hs.subset_extend (Set.subset_univ s) h2mem⟩
  let vals : hs.extend (Set.subset_univ s) -> ℝ :=
    fun i => if (i : ℝ) = Real.sqrt 2 then -Real.sqrt 2 else (i : ℝ)
  let L : ℝ →ₗ[ℚ] ℝ := B.constr ℚ vals
  have hB1 : B i1 = (1 : ℝ) := by
    simpa [B, i1] using Module.Basis.extend_apply_self hs i1
  have hB2 : B i2 = Real.sqrt 2 := by
    simpa [B, i2] using Module.Basis.extend_apply_self hs i2
  have hL1 : L 1 = 1 := by
    calc
      L 1 = L (B i1) := by rw [hB1]
      _ = vals i1 := by exact Module.Basis.constr_basis (b := B) ℚ vals i1
      _ = 1 := by
        have hne : ((i1 : hs.extend (Set.subset_univ s)) : ℝ) ≠ Real.sqrt 2 := by
          intro h
          have : (1 : ℝ) = Real.sqrt 2 := by simpa [i1] using h
          exact (ne_of_lt Real.one_lt_sqrt_two) this
        simp [vals, hne, i1]
  have hL2 : L (Real.sqrt 2) = -Real.sqrt 2 := by
    calc
      L (Real.sqrt 2) = L (B i2) := by rw [hB2]
      _ = vals i2 := by exact Module.Basis.constr_basis (b := B) ℚ vals i2
      _ = -Real.sqrt 2 := by simp [vals, i2]
  refine ⟨L, ?_, ?_, ?_⟩
  · intro x y
    exact L.map_add x y
  · have hm := L.map_smul ((2 : ℚ) / 3) (1 : ℝ)
    rw [hL1] at hm
    norm_num [Rat.smul_def] at hm ⊢
    exact hm
  · refine ⟨Real.sqrt 2, 0, ?_, ?_⟩
    · exact Real.sqrt_pos.2 (by norm_num)
    · rw [hL2]
      simp

/-- Applying an additive real pathology after centering. -/
noncomputable def badRankCentered3 (phi : ℝ -> ℝ)
    (x : Fin 3 -> ℝ) : Fin 3 -> ℝ :=
  fun i => phi (center3 x i)

/-- Graph label: centering an additive real pathology gives perturbation additivity. -/
lemma badRankCentered3_additive {phi : ℝ -> ℝ}
    (hadd : ∀ x y : ℝ, phi (x + y) = phi x + phi y) :
    PerturbationAdditive (badRankCentered3 phi) := by
  intro x y
  ext i
  have hcenter :
      center3 (x + y) i = center3 x i + center3 y i := by
    simp [center3]
    rw [Finset.sum_add_distrib]
    ring_nf
  simp [badRankCentered3, hcenter, hadd]

/-- Graph label: centered real pathologies are relabeling equivariant. -/
lemma badRankCentered3_equivariant {phi : ℝ -> ℝ} :
    RelabelingEquivariant (badRankCentered3 phi) := by
  intro sigma x
  ext i
  have hsum :
      (∑ j : Fin 3, x (sigma j)) = ∑ j : Fin 3, x j := by
    exact Equiv.sum_comp sigma x
  simp [badRankCentered3, center3, permute, hsum]

/-- Graph label: centered real pathologies are fully scale invariant. -/
lemma badRankCentered3_scale_invariant_full {phi : ℝ -> ℝ} :
    ScaleInvariantFull (badRankCentered3 phi) := by
  intro x c
  ext i
  have hsum :
      (∑ j : Fin 3, (x j + c)) = (∑ j : Fin 3, x j) + 3 * c := by
    norm_num [Fin.sum_univ_three]
    ring
  simp [badRankCentered3, center3, hsum]
  ring_nf

/-- Graph label: a calibrated real pathology gives a calibrated centered transform. -/
lemma badRankCentered3_calibrated {phi : ℝ -> ℝ}
    (hcal : phi ((2 : ℝ) / 3) = (2 : ℝ) / 3) :
    NaturallyCalibrated 1 (badRankCentered3 phi) := by
  have h20 : (2 : Fin 3) ≠ 0 := by decide
  norm_num [NaturallyCalibrated, badRankCentered3, center3, basis, first,
    Fin.sum_univ_three, h20, hcal]

/-- Graph label: a non-order-preserving real pathology violates rank monotonicity. -/
lemma badRankCentered3_not_rank_monotone {phi : ℝ -> ℝ}
    (hbad : ∃ a b : ℝ, a > b ∧ ¬ phi a > phi b) :
    ¬ RankMonotone (badRankCentered3 phi) := by
  rintro hmono
  obtain ⟨a, b, hab, hbad'⟩ := hbad
  let x : Fin 3 -> ℝ :=
    fun i => if i = (⟨0, by norm_num⟩ : Fin 3) then a
      else if i = (⟨1, by norm_num⟩ : Fin 3) then b
      else -(a + b)
  have hx : x (⟨0, by norm_num⟩ : Fin 3) >
      x (⟨1, by norm_num⟩ : Fin 3) := by
    simp [x, hab]
  have hT := (hmono x (⟨0, by norm_num⟩ : Fin 3)
    (⟨1, by norm_num⟩ : Fin 3)).mp hx
  have h10 : (1 : Fin 3) ≠ 0 := by decide
  have h20 : (2 : Fin 3) ≠ 0 := by decide
  have h21 : (2 : Fin 3) ≠ 1 := by decide
  norm_num [badRankCentered3, center3, x, Fin.sum_univ_three, h10, h20, h21] at hT
  exact hbad' hT

/-- Graph label: independence witness when rank monotonicity is omitted. -/
theorem rank_monotonicity_axiom_necessary :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      PerturbationAdditive T ∧ RelabelingEquivariant T ∧ ScaleInvariantFull T ∧
        NaturallyCalibrated 1 T ∧ ¬ RankMonotone T := by
  obtain ⟨phi, hadd, hcal, hbad⟩ := exists_rank_bad_additive_calibrated_phi
  exact ⟨badRankCentered3 phi, badRankCentered3_additive hadd,
    badRankCentered3_equivariant, badRankCentered3_scale_invariant_full,
    badRankCentered3_calibrated hcal, badRankCentered3_not_rank_monotone hbad⟩

/-- The five independence witnesses for the axioms in dimension 3. -/
structure AxiomNecessityWitnesses : Prop where
  withoutCalibration :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ RelabelingEquivariant T ∧
        ScaleInvariant T ∧ DiffersFromCLR T
  withoutScaleInvariance :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ RelabelingEquivariant T ∧
        NaturallyCalibrated 1 T ∧ DiffersFromCLR T
  withoutRelabelingEquivariance :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ PerturbationAdditive T ∧ ScaleInvariant T ∧
        NaturallyCalibrated 1 T ∧ DiffersFromCLR T
  withoutPerturbationAdditivity :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      RankMonotone T ∧ RelabelingEquivariant T ∧ ScaleInvariantFull T ∧
        NaturallyCalibrated 1 T ∧ ¬ PerturbationAdditive T
  withoutRankMonotonicity :
    ∃ T : (Fin 3 -> ℝ) -> (Fin 3 -> ℝ),
      PerturbationAdditive T ∧ RelabelingEquivariant T ∧ ScaleInvariantFull T ∧
        NaturallyCalibrated 1 T ∧ ¬ RankMonotone T

/--
Graph label: LaTeX tether for necessity of the axioms.

Named Lean statement for Supplementary Note Proposition `prop:axiom-necessity`.
It packages the five dimension-3 independence witnesses.
-/
theorem axiom_necessity_tether : AxiomNecessityWitnesses where
  withoutCalibration := calibration_axiom_necessary
  withoutScaleInvariance := scale_invariance_axiom_necessary
  withoutRelabelingEquivariance := relabeling_equivariance_axiom_necessary
  withoutPerturbationAdditivity := perturbation_additivity_axiom_necessary
  withoutRankMonotonicity := rank_monotonicity_axiom_necessary

end Invariance

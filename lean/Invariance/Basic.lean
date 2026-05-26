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
Additive characterization of CLR on the log-domain.

The formal dimension is `D + 2`, which is Lean's constructive way of saying
that the number of coordinates is at least two.  The result is the first
section of the supplement after passing to logarithmic coordinates: a
continuous additive, permutation-equivariant, scale-invariant, naturally
calibrated map is the centered log-ratio projection.
-/
theorem clr_characterization (D : Nat)
    (T : (Fin (D + 2) -> ℝ) -> (Fin (D + 2) -> ℝ))
    (hcont : Continuous T)
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

  let f : (Fin n -> ℝ) →+ (Fin n -> ℝ) := {
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
  let L : (Fin n -> ℝ) →L[ℝ] (Fin n -> ℝ) := f.toRealLinearMap hcont

  have hLT : ∀ x : Fin n -> ℝ, L x = T x := by
    intro x
    rfl
  have hpermL :
      ∀ sigma : Equiv.Perm (Fin n), ∀ x : Fin n -> ℝ,
        L (permute sigma x) = permute sigma (L x) := by
    intro sigma x
    simpa [L, f] using hperm sigma x

  let a : ℝ := L (basis i0) i0
  let b : ℝ := L (basis i0) i1

  have hEi : ∀ i : Fin n,
      L (basis i) = permute (Equiv.swap i0 i) (L (basis i0)) := by
    intro i
    have h := hpermL (Equiv.swap i0 i) (basis i0)
    simpa [permute_swap_basis] using h

  have hdiag : ∀ i : Fin n, L (basis i) i = a := by
    intro i
    have h := congrArg (fun v => v i) (hEi i)
    by_cases hi : i = i0
    · subst hi
      simp [a, permute] at h ⊢
    · simpa [a, permute, hi] using h

  have hOff0 : ∀ j : Fin n, j ≠ i0 -> L (basis i0) j = b := by
    intro j hj
    by_cases hj1 : j = i1
    · subst hj1
      rfl
    · have hfix :
          permute (Equiv.swap i1 j) (basis i0) = basis i0 := by
        exact permute_swap_basis_of_ne hi01 hj.symm
      have hper := hpermL (Equiv.swap i1 j) (basis i0)
      have hper0 : L (basis i0) =
          permute (Equiv.swap i1 j) (L (basis i0)) := by
        simpa [hfix] using hper
      have hcoord := congrArg (fun v => v i1) hper0
      have hb' : b = L (basis i0) j := by
        simpa [b, permute, hj1] using hcoord
      exact hb'.symm

  have hOff : ∀ i j : Fin n, i ≠ j -> L (basis i) j = b := by
    intro i j hij
    have h := congrArg (fun v => v j) (hEi i)
    by_cases hi : i = i0
    · subst hi
      simpa [permute] using hOff0 j hij.symm
    · have hne : Equiv.swap i0 i j ≠ i0 := by
        intro hs
        have h' := congrArg (Equiv.swap i0 i) hs
        have hswap0 : Equiv.swap i0 i i0 = i := by
          simp
        have hji : j = i := by
          simpa [hswap0] using h'
        exact hij hji.symm
      calc
        L (basis i) j = L (basis i0) (Equiv.swap i0 i j) := by
          simpa [permute] using h
        _ = b := hOff0 (Equiv.swap i0 i j) hne

  have hscaleeq : a + (D + 1 : ℝ) * b = 0 := by
    have hLone : L (1 : Fin n -> ℝ) = 0 := by
      simpa [L, f, n] using hscale (1 : ℝ)
    have hmap0 :
        L (∑ j : Fin n, basis j) i0 =
          (∑ j : Fin n, L (basis j) i0) := by
      simpa using congrArg (fun v => v i0)
        (map_sum L (fun j : Fin n => basis j) (Finset.univ : Finset (Fin n)))
    have hsumL_zero : (∑ j : Fin n, L (basis j) i0) = 0 := by
      have hleft : L (∑ j : Fin n, basis j) i0 = 0 := by
        simpa [sum_basis_eq_one] using congrArg (fun v => v i0) hLone
      linarith [hmap0, hleft]
    have hcoeff :
        (∑ j : Fin n, L (basis j) i0) =
          (∑ j : Fin n, (if j = i0 then a else b)) := by
      refine Finset.sum_congr rfl ?_
      intro j _
      by_cases hji : j = i0
      · subst hji
        simp [a, hdiag]
      · simp [hji, hOff j i0 hji]
    have hones : (∑ j : Fin n, (1 : ℝ)) = (D + 2 : ℝ) := by
      simp [n]
    have hsum_if :
        (∑ j : Fin n, (if j = i0 then a else b)) =
          a + (D + 1 : ℝ) * b := by
      calc
        (∑ j : Fin n, (if j = i0 then a else b))
            = (∑ j : Fin n, (1 : ℝ) * (if j = i0 then a else b)) := by
                simp
        _ = (a - b) * (1 : ℝ) + b * (∑ j : Fin n, (1 : ℝ)) := by
              simpa using sum_mul_if_eq_add (fun _ : Fin n => (1 : ℝ)) i0 a b
        _ = a + (D + 1 : ℝ) * b := by
              rw [hones]
              ring
    linarith [hsumL_zero, hcoeff, hsum_if]

  have hcalL : a = (D + 1 : ℝ) / (D + 2 : ℝ) := by
    simpa [a, i0, L, f, n] using hcal

  have hb : b = -(1 / (D + 2 : ℝ)) := by
    have hmain : (D + 1 : ℝ) * (1 / (D + 2 : ℝ) + b) = 0 := by
      calc
        (D + 1 : ℝ) * (1 / (D + 2 : ℝ) + b)
            = (D + 1 : ℝ) / (D + 2 : ℝ) + (D + 1 : ℝ) * b := by
                ring
        _ = a + (D + 1 : ℝ) * b := by
              simp [hcalL]
        _ = 0 := hscaleeq
    have hsum : 1 / (D + 2 : ℝ) + b = 0 := by
      exact (mul_eq_zero.mp hmain).resolve_left hdim1_ne
    linarith

  have halpha : a - b = 1 := by
    rw [hcalL, hb]
    field_simp [hdim_ne]
    ring

  intro x
  ext i
  calc
    T x i = L x i := by rfl
    _ = L (∑ j : Fin n, x j • basis j) i := by
          rw [sum_smul_basis]
    _ = (∑ j : Fin n, L (x j • basis j) i) := by
          simpa using congrArg (fun v => v i)
            (map_sum L (fun j : Fin n => x j • basis j)
              (Finset.univ : Finset (Fin n)))
    _ = (∑ j : Fin n, x j * L (basis j) i) := by
          simp
    _ = (∑ j : Fin n, x j * (if j = i then a else b)) := by
          refine Finset.sum_congr rfl ?_
          intro j _
          by_cases hji : j = i
          · subst hji
            simp [hdiag, a]
          · simp [hji, hOff j i hji]
    _ = (a - b) * x i + b * (∑ j : Fin n, x j) :=
          sum_mul_if_eq_add x i a b
    _ = x i - (∑ j : Fin n, x j) / (D + 2 : ℝ) := by
          rw [halpha, hb]
          ring
    _ = clr x i := by
          simp [clr, n]

end Invariance

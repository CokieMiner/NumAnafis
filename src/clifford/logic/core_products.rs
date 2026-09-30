//! Algebraic products and grade operations for [`CliffordNumber`].
//!
//! Blade arithmetic over active blades only: geometric, outer, and inner
//! products with Cayley-table acceleration, plus grade projections,
//! involutions, and the squared norm.

use alloc::vec::Vec;

use crate::types::Real;

use super::{BladeAlgebra, CliffordNumber};
impl CliffordNumber {
    // Products

    /// Geometric product `self × other`.
    #[must_use]
    pub fn geometric_mul(&self, other: &Self) -> Self {
        if self.gens == other.gens {
            self.geo_same(other)
        } else {
            let (u, pa, pb) = self.gens.union_with(&other.gens);
            self.reembed(u.clone(), &pa)
                .geo_same(&other.reembed(u, &pb))
        }
    }

    fn geo_same(&self, other: &Self) -> Self {
        debug_assert_eq!(
            self.gens, other.gens,
            "geo_same requires identical generator sets"
        );
        let limit = self.blade_count();
        let mut out = Self::zero_unchecked(self.gens.clone());

        let use_cache = !self.gens.cayley_signs.is_empty();
        let stride = 1 << self.gens.len();

        if limit <= 32 {
            let lhs_mask = self
                .nonzero_blades()
                .fold(0_u32, |bits, (blade, _)| bits | (1_u32 << blade));
            let rhs_mask = other
                .nonzero_blades()
                .fold(0_u32, |bits, (blade, _)| bits | (1_u32 << blade));
            let mut remaining_a = lhs_mask;
            while remaining_a != 0 {
                let a =
                    usize::try_from(remaining_a.trailing_zeros()).expect("inline blade fits usize");
                remaining_a &= remaining_a - 1;
                let mut remaining_b = rhs_mask;
                while remaining_b != 0 {
                    let b = usize::try_from(remaining_b.trailing_zeros())
                        .expect("inline blade fits usize");
                    remaining_b &= remaining_b - 1;
                    let factor = self.gens.cayley_signs[a * stride + b];
                    if factor == 0 {
                        continue;
                    }
                    let raw = self.coeff(a) * other.coeff(b);
                    let term = if factor < 0 { -raw } else { raw };
                    let target = &mut out.coeffs_mut_slice()[a ^ b];
                    *target += &term;
                }
            }
            return out;
        }
        let lhs_active: Vec<usize> = self
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();
        let rhs_active: Vec<usize> = other
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();

        for &a in &lhs_active {
            let lhs = &self.coeffs_slice()[a];
            for &b in &rhs_active {
                let rhs = &other.coeffs_slice()[b];
                let factor = if use_cache {
                    self.gens.cayley_signs[a * stride + b]
                } else {
                    BladeAlgebra::mul_blades(&self.gens, a, b).map_or(0, |(_, f)| f)
                };
                if factor == 0 {
                    continue;
                }
                let mask = a ^ b;
                let raw = lhs * rhs;
                let term = if factor < 0 { -raw } else { raw };
                let target = &mut out.coeffs_mut_slice()[mask];
                *target += &term;
            }
        }
        out
    }

    /// Outer (wedge) product `self ∧ other`.
    #[must_use]
    pub fn outer_product(&self, other: &Self) -> Self {
        if self.gens == other.gens {
            self.outer_same(other)
        } else {
            let (u, pa, pb) = self.gens.union_with(&other.gens);
            self.reembed(u.clone(), &pa)
                .outer_same(&other.reembed(u, &pb))
        }
    }

    fn outer_same(&self, other: &Self) -> Self {
        debug_assert_eq!(
            self.gens, other.gens,
            "outer_same requires identical generator sets"
        );
        let n = self.gens.len();
        let limit = self.blade_count();
        let mut out = Self::zero_unchecked(self.gens.clone());

        let use_cache = !self.gens.cayley_signs.is_empty();
        let stride = 1 << n;

        let lhs_active: Vec<usize> = self
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();
        let rhs_active: Vec<usize> = other
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();

        for &a in &lhs_active {
            let lhs = &self.coeffs_slice()[a];
            for &b in &rhs_active {
                if a & b != 0 {
                    continue;
                }
                let rhs = &other.coeffs_slice()[b];
                // Cayley cache sign is valid for the outer product when a & b == 0:
                // no shared generators ⇒ no metric contractions, so the geometric
                // product sign equals pure swap parity.
                let factor = if use_cache {
                    self.gens.cayley_signs[a * stride + b]
                } else {
                    let mut swaps = 0_u32;
                    for i in 0..n {
                        if (b >> i) & 1 == 1 {
                            swaps += (a >> (i + 1)).count_ones();
                        }
                    }
                    if swaps & 1 == 1 { -1 } else { 1 }
                };
                if factor == 0 {
                    continue;
                }
                let raw = lhs * rhs;
                let term = if factor < 0 { -raw } else { raw };
                let mask = a | b; // equivalent to a ^ b since overlap is 0
                let target = &mut out.coeffs_mut_slice()[mask];
                *target += &term;
            }
        }
        out
    }

    /// Left contraction: grade-`r` and grade-`s` terms contribute the
    /// grade-`s-r` part of their geometric product when `r <= s`.
    #[must_use]
    pub fn left_contraction(&self, other: &Self) -> Self {
        if self.gens == other.gens {
            self.inner_same(other)
        } else {
            let (u, pa, pb) = self.gens.union_with(&other.gens);
            self.reembed(u.clone(), &pa)
                .inner_same(&other.reembed(u, &pb))
        }
    }

    fn inner_same(&self, other: &Self) -> Self {
        debug_assert_eq!(
            self.gens, other.gens,
            "inner_same requires identical generator sets"
        );
        let limit = self.blade_count();
        let mut out = Self::zero_unchecked(self.gens.clone());

        let use_cache = !self.gens.cayley_signs.is_empty();
        let stride = 1 << self.gens.len();

        let lhs_active: Vec<usize> = self
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();
        let rhs_active: Vec<usize> = other
            .coeffs_slice()
            .iter()
            .enumerate()
            .take(limit)
            .filter(|item| !item.1.is_zero())
            .map(|item| item.0)
            .collect();

        for &a in &lhs_active {
            let lhs = &self.coeffs_slice()[a];
            let ga = a.count_ones();
            for &b in &rhs_active {
                let gb = b.count_ones();
                if gb < ga || (a & b) != a {
                    continue;
                }
                let rhs = &other.coeffs_slice()[b];
                let factor = if use_cache {
                    self.gens.cayley_signs[a * stride + b]
                } else {
                    BladeAlgebra::mul_blades(&self.gens, a, b).map_or(0, |(_, f)| f)
                };
                if factor == 0 {
                    continue;
                }
                let mask = a ^ b;
                let raw = lhs * rhs;
                let term = if factor < 0 { -raw } else { raw };
                let target = &mut out.coeffs_mut_slice()[mask];
                *target += &term;
            }
        }
        out
    }

    /// Grade-0 part of the geometric product `⟨AB⟩₀`.
    #[must_use]
    pub fn scalar_product(&self, other: &Self) -> Real {
        if self.gens != other.gens {
            return self.geometric_mul(other).coeff(0).clone();
        }
        let use_cache = !self.gens.cayley_signs.is_empty();
        let stride = 1 << self.gens.len();
        let mut result = Real::zero();
        for a in 0..self.blade_count() {
            let lhs = &self.coeffs_slice()[a];
            if lhs.is_zero() {
                continue;
            }
            let rhs = &other.coeffs_slice()[a];
            if rhs.is_zero() {
                continue;
            }
            let factor = if use_cache {
                self.gens.cayley_signs[a * stride + a]
            } else {
                BladeAlgebra::mul_blades(&self.gens, a, a).map_or(0, |(_, f)| f)
            };
            if factor == 0 {
                continue;
            }
            let term = lhs * rhs;
            if factor < 0 {
                result = &result - &term;
            } else {
                result = &result + &term;
            }
        }
        result
    }
    // Grade operations

    /// Project onto the grade-`k` subspace.
    #[must_use]
    pub fn grade(&self, k: u32) -> Self {
        let mut out = Self::zero_unchecked(self.gens.clone());
        for blade in 0..self.blade_count() {
            if blade.count_ones() == k {
                out.coeffs_mut_slice()[blade] = self.coeffs_slice()[blade].clone();
            }
        }
        out
    }

    /// Reverse `Ã`: multiply grade-g blades by `(−1)^(g(g−1)/2)`.
    #[must_use]
    pub fn reverse(&self) -> Self {
        let mut out = self.clone();
        for blade in 0..self.blade_count() {
            if blade.count_ones() % 4 >= 2 {
                let c = -&out.coeffs_slice()[blade];
                out.coeffs_mut_slice()[blade] = c;
            }
        }
        out
    }

    /// Grade involution `Â`: multiply grade-g blades by `(−1)^g`.
    #[must_use]
    pub fn grade_involution(&self) -> Self {
        let mut out = self.clone();
        for blade in 0..self.blade_count() {
            if blade.count_ones() % 2 == 1 {
                let c = -&out.coeffs_slice()[blade];
                out.coeffs_mut_slice()[blade] = c;
            }
        }
        out
    }

    /// Clifford conjugate `A†`: reverse followed by grade involution.
    #[must_use]
    pub fn clifford_conjugate(&self) -> Self {
        let mut out = self.clone();
        for blade in 0..self.blade_count() {
            let g = blade.count_ones() % 4;
            let coefficient = self.coeffs_slice()[blade].clone();
            out.coeffs_mut_slice()[blade] = if g == 1 || g == 2 {
                -coefficient
            } else {
                coefficient
            };
        }
        out
    }

    /// Conjugate quadratic form `⟨A A†⟩₀`, where `A†` is the Clifford conjugate.
    ///
    /// On the CGA catalog subalgebras this is `a² + b²` for `a + b ci`,
    /// the sum of four squares for quaternions, `a² − b²` for `a + b sj`,
    /// and `a²` for `a + b eps`. For general multivectors it is a signed
    /// scalar quadratic form, not a positive distance or a general inverse.
    /// The reversion product is `self.scalar_product(&self.reverse())`.
    #[must_use]
    pub fn norm_sq(&self) -> Real {
        self.scalar_product(&self.clifford_conjugate())
    }
}

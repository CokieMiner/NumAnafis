use alloc::vec;

use crate::types::{RationalType, Real};

use super::{CliffordNumber, E_MINUS, E_PLUS, E1, E2, E3, GeneratorSet};

/// Conformal Geometric Algebra Cl(4,1): canonical generators, derived
/// elements (conformal points, quaternions, pseudoscalars), and aliases.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Cga;

impl Cga {
    // CGA Cl(4,1) full generator set

    /// The full CGA generator set: `(e1, e2, e3, e+, e−)` with signature `(4,1)`.
    ///
    /// This is the canonical generator set for Conformal Geometric Algebra.
    /// All derived constructions (quaternions, conformal points, etc.) live in
    /// this algebra.
    #[must_use]
    pub fn gens() -> GeneratorSet {
        GeneratorSet::from_validated(vec![
            (E1, 1_i8),
            (E2, 1_i8),
            (E3, 1_i8),
            (E_PLUS, 1_i8),
            (E_MINUS, -1_i8),
        ])
    }
    // Base generators — each is a grade-1 element of the full CGA algebra

    /// Helper: build a CGA multivector with a single generator blade set to 1.
    fn cga_unit(bit: usize) -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        mv.coeffs_mut_slice()[1_usize << bit] = Real::one();
        mv
    }

    /// Euclidean basis vector `e1` with `e1² = +1`.
    #[must_use]
    #[inline]
    pub fn e1() -> CliffordNumber {
        Self::cga_unit(0) // bit 0 → blade mask 0b00001
    }

    /// Euclidean basis vector `e2` with `e2² = +1`.
    #[must_use]
    #[inline]
    pub fn e2() -> CliffordNumber {
        Self::cga_unit(1) // bit 1 → blade mask 0b00010
    }

    /// Euclidean basis vector `e3` with `e3² = +1`.
    #[must_use]
    #[inline]
    pub fn e3() -> CliffordNumber {
        Self::cga_unit(2) // bit 2 → blade mask 0b00100
    }

    /// Extra positive dimension `e+` with `e+² = +1`.
    #[must_use]
    #[inline]
    pub fn e_plus() -> CliffordNumber {
        Self::cga_unit(3) // bit 3 → blade mask 0b01000
    }

    /// Extra negative dimension `e−` with `e−² = −1`.
    #[must_use]
    #[inline]
    pub fn e_minus() -> CliffordNumber {
        Self::cga_unit(4) // bit 4 → blade mask 0b10000
    }
    // Derived elements of the full CGA algebra.

    /// Conformal origin point: `eₒ = ½(e− − e+)`.
    ///
    /// Properties: `eₒ² = 0`, `eₒ · e∞ = −1`, and `eₒ e∞ + e∞ eₒ = −2`.
    #[must_use]
    #[inline]
    pub fn orig() -> CliffordNumber {
        (Self::e_minus() - Self::e_plus()) * Real::from_rational(RationalType::from_parts(1, 2))
    }

    /// Conformal infinity: `e∞ = e− + e+`.
    ///
    /// Properties: `e∞² = 0` and `eₒ · e∞ = −1`.
    /// Its span with the scalar identity is a dual-number subalgebra.
    #[must_use]
    #[inline]
    pub fn inf() -> CliffordNumber {
        Self::e_minus() + Self::e_plus()
    }

    /// Quaternion imaginary unit `i = −e2·e3` (grade-2 bivector).
    ///
    /// Satisfies Hamilton's conventions: `i² = j² = k² = ijk = −1`,
    /// `ij = k`, `jk = i`, and `ki = j`; reversing a pair changes the sign.
    #[must_use]
    #[inline]
    pub fn qi() -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        // e2∧e3 → bits 1,2 → blade mask 0b00110 = 6
        mv.coeffs_mut_slice()[6] = -Real::one();
        mv
    }

    /// Quaternion imaginary unit `j = −e3·e1` (grade-2 bivector).
    ///
    /// Together with [`Self::qi`] and [`Self::qk`], satisfies `ij = k`,
    /// `jk = i`, `ki = j`, and `i² = j² = k² = ijk = −1`.
    #[must_use]
    #[inline]
    pub fn qj() -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        // e3∧e1 → bits 0,2 → blade mask 0b00101 = 5
        mv.coeffs_mut_slice()[5] = Real::one();
        mv
    }

    /// Quaternion imaginary unit `k = −e1·e2` (grade-2 bivector).
    ///
    /// Together with [`Self::qi`] and [`Self::qj`], satisfies `ij = k`,
    /// `jk = i`, `ki = j`, and `i² = j² = k² = ijk = −1`.
    #[must_use]
    #[inline]
    pub fn qk() -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        // e1∧e2 → bits 0,1 → blade mask 0b00011 = 3
        mv.coeffs_mut_slice()[3] = -Real::one();
        mv
    }

    /// 3D pseudoscalar `I₃ = e1·e2·e3` (grade-3 trivector).
    ///
    /// Properties: `I₃² = −1`; commutes with all elements of the Euclidean
    /// sub-algebra, but anticommutes with `e+` and `e−`.
    /// It is the oriented Euclidean volume element; `−I₃ e1`, `−I₃ e2`,
    /// and `−I₃ e3` are the three quaternion units.
    #[must_use]
    #[inline]
    pub fn pseudo3d() -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        // e1∧e2∧e3 → bits 0,1,2 → blade mask 0b00111 = 7
        mv.coeffs_mut_slice()[7] = Real::one();
        mv
    }

    /// 5D pseudoscalar `I₅ = e1·e2·e3·e+·e−` (grade-5 blade).
    ///
    /// Properties: `I₅² = −1` (signature dependent: (−1)^{10} from grade
    /// times the metric product `(+1)⁴·(−1)`).
    /// It commutes with every CGA blade, so `a + b I₅` embeds ordinary
    /// complex numbers centrally into this real Clifford algebra.
    #[must_use]
    #[inline]
    pub fn pseudo5d() -> CliffordNumber {
        let gens = Self::gens();
        let mut mv = CliffordNumber::zero_unchecked(gens);
        // all 5 bits → blade mask 0b11111 = 31
        mv.coeffs_mut_slice()[31] = Real::one();
        mv
    }
    // Aliases — convenient names for derived CGA elements

    /// Dual / nilpotent unit `ε` — alias for [`Self::inf`] (`e∞ = e− + e+`).
    ///
    /// `ε² = 0` and `ε ≠ 0`: `a + b ε` obeys dual-number arithmetic.
    /// It commutes with Euclidean bivectors and anticommutes with Euclidean
    /// vectors. It is not a central dual unit for the whole CGA.
    #[must_use]
    #[inline]
    pub fn eps() -> CliffordNumber {
        Self::inf()
    }

    /// Central complex imaginary unit `i` — alias for [`Self::pseudo5d`].
    ///
    /// `ci² = −1` and `ci` commutes with every CGA blade, preserving the
    /// multiplication laws of ordinary complex scalar numbers.
    #[must_use]
    #[inline]
    pub fn ci() -> CliffordNumber {
        Self::pseudo5d()
    }

    /// Split-complex unit `j` — alias for [`Self::e_plus`] (`e+`).
    ///
    /// `j² = +1`: `a + b j` obeys split-complex (hyperbolic) arithmetic,
    /// with conjugate norm `a² − b²` and zero divisors `1 ± j`.
    /// It commutes with Euclidean bivectors and anticommutes with Euclidean
    /// vectors; cross-subalgebra products follow the CGA geometric product.
    #[must_use]
    #[inline]
    pub fn sj() -> CliffordNumber {
        Self::e_plus()
    }
}

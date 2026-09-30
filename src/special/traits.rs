//! Backend abstraction traits for special-function evaluation.
//!
//! [`SpecInt`] abstracts the integer orders (`n`, `l`, `m`) and [`SpecFloat`]
//! abstracts IEEE float operations, range-reduction constants, and the
//! precision-specific coefficient tables behind compilation walls.

use core::{
    fmt::Debug,
    ops::{Add, Div, Mul, Neg, Rem, Sub},
};

// ============================================================================
// SpecInt — integer abstraction
// ============================================================================

/// Trait abstracting integers for special function polynomial evaluations.
pub trait SpecInt:
    Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Debug
    + 'static
{
    fn zero() -> Self;
    fn one() -> Self;
    #[must_use]
    fn abs(self) -> Self;
    fn is_zero(self) -> bool;
    fn is_negative(self) -> bool;
    fn to_usize(self) -> usize;
    fn from_usize(v: usize) -> Self;
}

// ============================================================================
// SpecFloat — float abstraction
// ============================================================================

/// Trait abstracting floating point for special function implementations.
pub trait SpecFloat:
    Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + 'static
{
    type Int: SpecInt;

    fn eps() -> Self;
    const ERF_TERMS: usize;
    const DIGAMMA_SHIFT: usize;
    const LAMBERT_MAX_ITERATIONS: usize;
    const ZETA_BORWEIN_N: usize;
    const ZETA_EM_TERMS: usize;

    fn zero() -> Self;
    fn one() -> Self;
    #[must_use]
    fn neg_one() -> Self {
        -Self::one()
    }
    #[must_use]
    fn two() -> Self {
        Self::one() + Self::one()
    }
    #[must_use]
    fn half() -> Self {
        Self::one() / Self::two()
    }
    fn pi() -> Self;
    fn e() -> Self;
    fn frac_2_pi() -> Self;
    fn infinity() -> Self;
    #[must_use]
    fn neg_infinity() -> Self {
        -Self::infinity()
    }
    fn nan() -> Self;
    fn max_value() -> Self;

    // Cody-Waite range reduction constants for π/4 and 3π/4.
    fn pio4_hi() -> Self;
    fn pio4_lo() -> Self;
    fn pio34_hi() -> Self;
    fn pio34_lo() -> Self;

    fn from_int<I: SpecInt>(v: I) -> Self {
        let a = v.abs();
        let f = Self::from_usize(a.to_usize());
        if v.is_negative() { -f } else { f }
    }
    fn from_usize(v: usize) -> Self;

    #[must_use]
    fn abs(self) -> Self;
    #[must_use]
    fn signum(self) -> Self;
    #[must_use]
    fn sqrt(self) -> Self;
    #[must_use]
    fn sin(self) -> Self;
    #[must_use]
    fn cos(self) -> Self;
    #[must_use]
    fn atan2(self, other: Self) -> Self;
    #[must_use]
    fn max(self, other: Self) -> Self;
    #[must_use]
    fn mul_add(self, a: Self, b: Self) -> Self;
    #[must_use]
    fn exp(self) -> Self;
    #[must_use]
    fn ln(self) -> Self;
    #[must_use]
    fn fract(self) -> Self;
    #[must_use]
    fn round(self) -> Self;
    #[must_use]
    fn floor(self) -> Self;
    #[must_use]
    fn powf(self, exp: Self) -> Self;
    #[must_use]
    fn pow_int<I: SpecInt>(self, n: I) -> Self;

    fn is_nan(self) -> bool;
    fn is_infinite(self) -> bool;
    fn is_sign_negative(self) -> bool;
    fn is_sign_positive(self) -> bool {
        !self.is_sign_negative() && !self.is_nan()
    }

    // --- Coefficient arrays for special functions ---
    // Each backend stores these as const arrays in native precision.

    fn lanczos_coeffs() -> &'static [Self];
    fn bernoulli_pairs() -> &'static [(Self, Self)];

    fn stieltjes_coeffs() -> &'static [Self];
    fn stirling_coeffs() -> &'static [Self];
    fn zeta_ints() -> &'static [Self];

    fn besselj0_num_coeffs() -> &'static [Self];
    fn besselj0_den_coeffs() -> &'static [Self];
    fn besselj0_pcos_coeffs() -> &'static [Self];
    fn besselj0_psin_coeffs() -> &'static [Self];

    fn besselj1_num_coeffs() -> &'static [Self];
    fn besselj1_den_coeffs() -> &'static [Self];
    fn besselj1_pcos_coeffs() -> &'static [Self];
    fn besselj1_psin_coeffs() -> &'static [Self];

    fn bessely0_num_coeffs() -> &'static [Self];
    fn bessely0_den_coeffs() -> &'static [Self];
    fn bessely0_pcos_coeffs() -> &'static [Self];
    fn bessely0_psin_coeffs() -> &'static [Self];

    fn bessely1_num_coeffs() -> &'static [Self];
    fn bessely1_den_coeffs() -> &'static [Self];
    fn bessely1_pcos_coeffs() -> &'static [Self];
    fn bessely1_psin_coeffs() -> &'static [Self];

    fn besseli0_small_coeffs() -> &'static [Self];
    fn besseli0_large_coeffs() -> &'static [Self];

    fn besseli1_small_coeffs() -> &'static [Self];
    fn besseli1_large_coeffs() -> &'static [Self];

    fn besselk0_small_coeffs() -> &'static [Self];
    fn besselk0_large_coeffs() -> &'static [Self];

    fn besselk1_small_coeffs() -> &'static [Self];
    fn besselk1_large_coeffs() -> &'static [Self];

    fn besselj_split() -> Self;
    fn bessel_miller_seed() -> Self;

    fn besselj0_root1() -> Self;
    fn besselj0_root2() -> Self;
    fn besselj1_root1() -> Self;
    fn besselj1_root2() -> Self;
}

// ============================================================================
// SpecFloat / SpecInt trait implementations for f32/i32 backend
// ============================================================================

#[cfg(not(feature = "std"))]
use libm::{
    atan2, atan2f, cos, cosf, exp, expf, floor, floorf, fma, fmaf, log, logf, pow, powf, round,
    roundf, sin, sinf, sqrt, sqrtf, trunc, truncf,
};

use super::{SpecFloat, SpecInt, coeffs_f32, coeffs_f64};

/// Implements the [`SpecInt`] trait for a concrete integer type.
macro_rules! impl_spec_int {
    ($ty:ty) => {
        impl SpecInt for $ty {
            #[inline]
            fn zero() -> Self {
                0
            }
            #[inline]
            fn one() -> Self {
                1
            }
            #[inline]
            fn abs(self) -> Self {
                <$ty>::abs(self)
            }
            #[inline]
            fn is_zero(self) -> bool {
                self == 0
            }
            #[inline]
            fn is_negative(self) -> bool {
                self < 0
            }
            #[inline]
            fn to_usize(self) -> usize {
                usize::try_from(self).unwrap_or(0)
            }
            #[inline]
            fn from_usize(v: usize) -> Self {
                <$ty>::try_from(v).unwrap_or(Self::MAX)
            }
        }
    };
}

/// Implements the [`SpecFloat`] trait for a concrete float type.
///
/// Precision-specific scalars and coefficient tables live in the sibling
/// `coeffs_f32` / `coeffs_f64` modules; the macro only wires `std`/*libm*
/// dispatch. Call sites stay small and diffable.
macro_rules! impl_spec_float {
    (
        $ty:ty,
        $int:ty,
        $coeffs:ident,
        sqrt_fn = $sqrt_fn:path,
        sin_fn = $sin_fn:path,
        cos_fn = $cos_fn:path,
        atan2_fn = $atan2_fn:path,
        exp_fn = $exp_fn:path,
        ln_fn = $ln_fn:path,
        trunc_fn = $trunc_fn:path,
        round_fn = $round_fn:path,
        floor_fn = $floor_fn:path,
        pow_fn = $pow_fn:path,
        mul_add_fn = $mul_add_fn:path,
    ) => {
        impl SpecFloat for $ty {
            type Int = $int;

            #[inline]
            fn eps() -> Self {
                <$ty>::EPSILON
            }
            const ERF_TERMS: usize = $coeffs::ERF_TERMS;
            const DIGAMMA_SHIFT: usize = $coeffs::DIGAMMA_SHIFT;
            const LAMBERT_MAX_ITERATIONS: usize = $coeffs::LAMBERT_MAX_ITERATIONS;
            const ZETA_BORWEIN_N: usize = $coeffs::ZETA_BORWEIN_N;
            const ZETA_EM_TERMS: usize = $coeffs::ZETA_EM_TERMS;

            #[inline]
            fn zero() -> Self {
                0.0
            }
            #[inline]
            fn one() -> Self {
                1.0
            }
            #[inline]
            fn pi() -> Self {
                $coeffs::PI
            }
            #[inline]
            fn e() -> Self {
                $coeffs::E
            }
            #[inline]
            fn frac_2_pi() -> Self {
                $coeffs::FRAC_2_PI
            }
            #[inline]
            fn infinity() -> Self {
                <$ty>::INFINITY
            }
            #[inline]
            fn nan() -> Self {
                <$ty>::NAN
            }
            #[inline]
            fn max_value() -> Self {
                <$ty>::MAX
            }

            // Cody-Waite range reduction constants for π/4 and 3π/4
            #[inline]
            fn pio4_hi() -> Self {
                $coeffs::PIO4_HI
            }
            #[inline]
            fn pio4_lo() -> Self {
                $coeffs::PIO4_LO
            }
            #[inline]
            fn pio34_hi() -> Self {
                $coeffs::PIO34_HI
            }
            #[inline]
            fn pio34_lo() -> Self {
                $coeffs::PIO34_LO
            }

            #[inline]
            fn from_usize(v: usize) -> Self {
                #[expect(clippy::as_conversions, reason = "Primitive usize to float conversion")]
                #[expect(
                    clippy::cast_precision_loss,
                    reason = "Unavoidable precision loss for huge integers"
                )]
                {
                    v as Self
                }
            }

            #[inline]
            fn abs(self) -> Self {
                self.abs()
            }
            #[inline]
            fn signum(self) -> Self {
                self.signum()
            }
            // Std-only inherent methods route to libm without std. The cfg is
            // explicit (rather than a trait import) because inside this impl,
            // a trait fallback would resolve to the method itself and recurse.
            #[inline]
            fn sqrt(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.sqrt()
                }
                #[cfg(not(feature = "std"))]
                {
                    $sqrt_fn(self)
                }
            }
            #[inline]
            fn sin(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.sin()
                }
                #[cfg(not(feature = "std"))]
                {
                    $sin_fn(self)
                }
            }
            #[inline]
            fn cos(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.cos()
                }
                #[cfg(not(feature = "std"))]
                {
                    $cos_fn(self)
                }
            }
            #[inline]
            fn atan2(self, other: Self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.atan2(other)
                }
                #[cfg(not(feature = "std"))]
                {
                    $atan2_fn(self, other)
                }
            }
            #[inline]
            fn max(self, other: Self) -> Self {
                <$ty>::max(self, other)
            }
            #[inline]
            fn mul_add(self, a: Self, b: Self) -> Self {
                #[cfg(feature = "std")]
                {
                    <$ty>::mul_add(self, a, b)
                }
                // Both are single-rounding fused multiply-add.
                #[cfg(not(feature = "std"))]
                {
                    $mul_add_fn(self, a, b)
                }
            }
            #[inline]
            fn exp(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.exp()
                }
                #[cfg(not(feature = "std"))]
                {
                    $exp_fn(self)
                }
            }
            #[inline]
            fn ln(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.ln()
                }
                #[cfg(not(feature = "std"))]
                {
                    $ln_fn(self)
                }
            }
            #[inline]
            fn fract(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.fract()
                }
                // libm provides no fract: inherent is self - self.trunc().
                #[cfg(not(feature = "std"))]
                {
                    self - $trunc_fn(self)
                }
            }
            #[inline]
            fn round(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.round()
                }
                #[cfg(not(feature = "std"))]
                {
                    $round_fn(self)
                }
            }
            #[inline]
            fn floor(self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.floor()
                }
                #[cfg(not(feature = "std"))]
                {
                    $floor_fn(self)
                }
            }
            #[inline]
            fn powf(self, exp: Self) -> Self {
                #[cfg(feature = "std")]
                {
                    self.powf(exp)
                }
                #[cfg(not(feature = "std"))]
                {
                    $pow_fn(self, exp)
                }
            }
            #[inline]
            fn pow_int<I: SpecInt>(self, n: I) -> Self {
                let mut result = Self::one();
                let mut base = self;
                let mut e = n.abs();
                let two = I::from_usize(2);
                while !e.is_zero() {
                    if (e % two) != I::zero() {
                        result *= base;
                    }
                    base *= base;
                    e = e / two;
                }
                if n.is_negative() {
                    Self::one() / result
                } else {
                    result
                }
            }
            #[inline]
            fn is_nan(self) -> bool {
                <$ty>::is_nan(self)
            }
            #[inline]
            fn is_infinite(self) -> bool {
                <$ty>::is_infinite(self)
            }
            #[inline]
            fn is_sign_negative(self) -> bool {
                <$ty>::is_sign_negative(self)
            }

            // --- Coefficient arrays — stored as native-precision consts ---

            #[inline]
            fn lanczos_coeffs() -> &'static [Self] {
                $coeffs::LANCZOS
            }
            #[inline]
            fn bernoulli_pairs() -> &'static [(Self, Self)] {
                $coeffs::BERNOULLI
            }
            #[inline]
            fn stieltjes_coeffs() -> &'static [Self] {
                $coeffs::STIELTJES
            }
            #[inline]
            fn stirling_coeffs() -> &'static [Self] {
                $coeffs::STIRLING
            }
            #[inline]
            fn zeta_ints() -> &'static [Self] {
                $coeffs::ZETA_INTS
            }

            #[inline]
            fn besselj0_num_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ0_NUM
            }
            #[inline]
            fn besselj0_den_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ0_DEN
            }
            #[inline]
            fn besselj0_pcos_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ0_PCOS
            }
            #[inline]
            fn besselj0_psin_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ0_PSIN
            }

            #[inline]
            fn besselj1_num_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ1_NUM
            }
            #[inline]
            fn besselj1_den_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ1_DEN
            }
            #[inline]
            fn besselj1_pcos_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ1_PCOS
            }
            #[inline]
            fn besselj1_psin_coeffs() -> &'static [Self] {
                $coeffs::BESSELJ1_PSIN
            }

            #[inline]
            fn bessely0_num_coeffs() -> &'static [Self] {
                $coeffs::BESSELY0_NUM
            }
            #[inline]
            fn bessely0_den_coeffs() -> &'static [Self] {
                $coeffs::BESSELY0_DEN
            }
            #[inline]
            fn bessely0_pcos_coeffs() -> &'static [Self] {
                $coeffs::BESSELY0_PCOS
            }
            #[inline]
            fn bessely0_psin_coeffs() -> &'static [Self] {
                $coeffs::BESSELY0_PSIN
            }

            #[inline]
            fn bessely1_num_coeffs() -> &'static [Self] {
                $coeffs::BESSELY1_NUM
            }
            #[inline]
            fn bessely1_den_coeffs() -> &'static [Self] {
                $coeffs::BESSELY1_DEN
            }
            #[inline]
            fn bessely1_pcos_coeffs() -> &'static [Self] {
                $coeffs::BESSELY1_PCOS
            }
            #[inline]
            fn bessely1_psin_coeffs() -> &'static [Self] {
                $coeffs::BESSELY1_PSIN
            }

            #[inline]
            fn besseli0_small_coeffs() -> &'static [Self] {
                $coeffs::BESSELI0_SMALL
            }
            #[inline]
            fn besseli0_large_coeffs() -> &'static [Self] {
                $coeffs::BESSELI0_LARGE
            }

            #[inline]
            fn besseli1_small_coeffs() -> &'static [Self] {
                $coeffs::BESSELI1_SMALL
            }
            #[inline]
            fn besseli1_large_coeffs() -> &'static [Self] {
                $coeffs::BESSELI1_LARGE
            }

            #[inline]
            fn besselk0_small_coeffs() -> &'static [Self] {
                $coeffs::BESSELK0_SMALL
            }
            #[inline]
            fn besselk0_large_coeffs() -> &'static [Self] {
                $coeffs::BESSELK0_LARGE
            }

            #[inline]
            fn besselk1_small_coeffs() -> &'static [Self] {
                $coeffs::BESSELK1_SMALL
            }
            #[inline]
            fn besselk1_large_coeffs() -> &'static [Self] {
                $coeffs::BESSELK1_LARGE
            }

            #[inline]
            fn besselj_split() -> Self {
                $coeffs::BESSELJ_SPLIT
            }
            #[inline]
            fn bessel_miller_seed() -> Self {
                $coeffs::BESSEL_MILLER_SEED
            }

            #[inline]
            fn besselj0_root1() -> Self {
                $coeffs::BESSELJ0_ROOT1
            }
            #[inline]
            fn besselj0_root2() -> Self {
                $coeffs::BESSELJ0_ROOT2
            }
            #[inline]
            fn besselj1_root1() -> Self {
                $coeffs::BESSELJ1_ROOT1
            }
            #[inline]
            fn besselj1_root2() -> Self {
                $coeffs::BESSELJ1_ROOT2
            }
        }
    };
}

impl_spec_int!(i32);

impl SpecInt for usize {
    #[inline]
    fn zero() -> Self {
        0
    }
    #[inline]
    fn one() -> Self {
        1
    }
    #[inline]
    fn abs(self) -> Self {
        self
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
    #[inline]
    fn is_negative(self) -> bool {
        false
    }
    #[inline]
    fn to_usize(self) -> usize {
        self
    }
    #[inline]
    fn from_usize(v: usize) -> Self {
        v
    }
}

impl_spec_float!(
    f32,
    i32,
    coeffs_f32,
    sqrt_fn = sqrtf,
    sin_fn = sinf,
    cos_fn = cosf,
    atan2_fn = atan2f,
    exp_fn = expf,
    ln_fn = logf,
    trunc_fn = truncf,
    round_fn = roundf,
    floor_fn = floorf,
    pow_fn = powf,
    mul_add_fn = fmaf,
);

// ============================================================================
// SpecFloat / SpecInt trait implementations for f64/i64 backend
// ============================================================================

impl_spec_int!(i64);

impl_spec_float!(
    f64,
    i64,
    coeffs_f64,
    sqrt_fn = sqrt,
    sin_fn = sin,
    cos_fn = cos,
    atan2_fn = atan2,
    exp_fn = exp,
    ln_fn = log,
    trunc_fn = trunc,
    round_fn = round,
    floor_fn = floor,
    pow_fn = pow,
    mul_add_fn = fma,
);

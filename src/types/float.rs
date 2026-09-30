//! Native floating-point arithmetic and elementary kernels.

#![expect(
    clippy::trivially_copy_pass_by_ref,
    reason = "Backend operands are borrowed so owned approximate representations can reuse storage"
)]

use core::{
    cmp::Ordering,
    f64::consts::{E, PI},
};

use libm::{ceil, floor, fma, hypot, modf, round, trunc};

#[cfg(all(not(feature = "std"), not(test)))]
use crate::ext::F64Ext;

use super::IntType;

// TODO(mp): supply owned float arithmetic and constants here; primitive projections
// remain explicit and must not be used by internal numerical kernels.
/// Native approximate representation.
pub type FloatType = f64;

/// Primitive arithmetic, elementary kernels, and machine constants for the float backend.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct FloatMath;

impl FloatMath {
    /// Explicit projection to the primitive float representation.
    pub const fn to_f64(value: &FloatType) -> f64 {
        *value
    }
    /// Backend NaN.
    pub const fn nan() -> FloatType {
        f64::NAN
    }
    /// Backend approximation to pi.
    pub const fn pi() -> FloatType {
        PI
    }
    /// Backend approximation to Euler's number.
    pub const fn e() -> FloatType {
        E
    }
    /// Total order of approximate backend values, including signed zeros and NaNs.
    pub fn total_cmp(lhs: &FloatType, rhs: &FloatType) -> Ordering {
        lhs.total_cmp(rhs)
    }
    /// Sign bit, including negative zero.
    pub const fn is_sign_negative(value: &FloatType) -> bool {
        value.is_sign_negative()
    }
    /// Approximate addition.
    pub fn add(lhs: &FloatType, rhs: &FloatType) -> FloatType {
        lhs + rhs
    }
    /// Approximate subtraction.
    pub fn sub(lhs: &FloatType, rhs: &FloatType) -> FloatType {
        lhs - rhs
    }
    /// Approximate multiplication.
    pub fn mul(lhs: &FloatType, rhs: &FloatType) -> FloatType {
        lhs * rhs
    }
    /// Approximate division.
    pub fn div(lhs: &FloatType, rhs: &FloatType) -> FloatType {
        lhs / rhs
    }
    /// Approximate negation.
    pub fn neg(value: &FloatType) -> FloatType {
        -value
    }
    /// Updates an approximate accumulator in place.
    pub fn add_assign(lhs: &mut FloatType, rhs: &FloatType) {
        *lhs += rhs;
    }
    /// Subtracts from an approximate accumulator in place.
    pub fn sub_assign(lhs: &mut FloatType, rhs: &FloatType) {
        *lhs -= rhs;
    }
    /// Scaled approximate complex division for finite nonzero denominators.
    /// Returns None when the algebraic exceptional-value path is required.
    pub fn complex_div(
        a: &FloatType,
        b: &FloatType,
        c: &FloatType,
        d: &FloatType,
    ) -> Option<(FloatType, FloatType)> {
        let scale = c.abs().max(d.abs());
        if !scale.is_finite() || scale == 0.0 {
            return None;
        }
        let (cs, ds) = (c / scale, d / scale);
        let denom = Self::mul_add(&cs, &cs, &(ds * ds));
        let (ascaled, bscaled) = (a / scale, b / scale);
        Some((
            Self::mul_add(&ascaled, &cs, &(bscaled * ds)) / denom,
            Self::mul_add(&bscaled, &cs, &-(ascaled * ds)) / denom,
        ))
    }

    /// Fused multiply-add in the active approximate backend.
    #[must_use]
    pub fn mul_add(a: &FloatType, b: &FloatType, c: &FloatType) -> FloatType {
        fma(*a, *b, *c)
    }

    /// Machine epsilon of the active backend.
    #[must_use]
    pub const fn epsilon() -> FloatType {
        f64::EPSILON
    }
    /// Convert an integer to the approximate backend, rounding when necessary.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "Explicit approximate projection rounds integers to f64"
    )]
    #[must_use]
    pub const fn from_int(value: &IntType) -> FloatType {
        *value as FloatType
    }
    /// Project a rational to the approximate backend.
    #[expect(
        clippy::as_conversions,
        clippy::cast_precision_loss,
        reason = "Explicit approximate projection rounds unsigned denominators to f64"
    )]
    #[must_use]
    pub fn from_ratio(numer: &IntType, denom: &u64) -> FloatType {
        Self::from_int(numer) / *denom as FloatType
    }
    /// Euclidean length, scaled to avoid squaring overflow.
    #[must_use]
    pub fn hypot(a: &FloatType, b: &FloatType) -> FloatType {
        hypot(*a, *b)
    }
    /// Backend sqrt kernel.
    #[must_use]
    pub fn sqrt(value: &FloatType) -> FloatType {
        value.sqrt()
    }
    /// Backend cbrt kernel.
    #[must_use]
    pub fn cbrt(value: &FloatType) -> FloatType {
        value.cbrt()
    }
    /// Backend exp kernel.
    #[must_use]
    pub fn exp(value: &FloatType) -> FloatType {
        value.exp()
    }
    /// Backend ln kernel.
    #[must_use]
    pub fn ln(value: &FloatType) -> FloatType {
        value.ln()
    }
    /// Backend sin kernel.
    #[must_use]
    pub fn sin(value: &FloatType) -> FloatType {
        value.sin()
    }
    /// Backend cos kernel.
    #[must_use]
    pub fn cos(value: &FloatType) -> FloatType {
        value.cos()
    }
    /// Backend tan kernel.
    #[must_use]
    pub fn tan(value: &FloatType) -> FloatType {
        value.tan()
    }
    /// Backend asin kernel.
    #[must_use]
    pub fn asin(value: &FloatType) -> FloatType {
        value.asin()
    }
    /// Backend acos kernel.
    #[must_use]
    pub fn acos(value: &FloatType) -> FloatType {
        value.acos()
    }
    /// Backend atan kernel.
    #[must_use]
    pub fn atan(value: &FloatType) -> FloatType {
        value.atan()
    }
    /// Backend sinh kernel.
    #[must_use]
    pub fn sinh(value: &FloatType) -> FloatType {
        value.sinh()
    }
    /// Backend cosh kernel.
    #[must_use]
    pub fn cosh(value: &FloatType) -> FloatType {
        value.cosh()
    }
    /// Backend tanh kernel.
    #[must_use]
    pub fn tanh(value: &FloatType) -> FloatType {
        value.tanh()
    }
    /// Backend asinh kernel.
    #[must_use]
    pub fn asinh(value: &FloatType) -> FloatType {
        value.asinh()
    }
    /// Backend acosh kernel.
    #[must_use]
    pub fn acosh(value: &FloatType) -> FloatType {
        value.acosh()
    }
    /// Backend atanh kernel.
    #[must_use]
    pub fn atanh(value: &FloatType) -> FloatType {
        value.atanh()
    }
    /// Backend floor kernel.
    #[must_use]
    pub fn floor(value: &FloatType) -> FloatType {
        floor(*value)
    }
    /// Backend ceil kernel.
    #[must_use]
    pub fn ceil(value: &FloatType) -> FloatType {
        ceil(*value)
    }
    /// Backend round kernel.
    #[must_use]
    pub fn round(value: &FloatType) -> FloatType {
        round(*value)
    }
    /// Backend trunc kernel.
    #[must_use]
    pub fn trunc(value: &FloatType) -> FloatType {
        trunc(*value)
    }
    /// Backend fract kernel.
    #[must_use]
    pub fn fract(value: &FloatType) -> FloatType {
        modf(*value).0
    }
    /// Backend abs kernel.
    #[must_use]
    pub const fn abs(value: &FloatType) -> FloatType {
        value.abs()
    }

    /// Backend expm1 kernel.
    #[must_use]
    pub fn expm1(value: &FloatType) -> FloatType {
        value.exp_m1()
    }
    /// Backend log1p kernel.
    #[must_use]
    pub fn log1p(value: &FloatType) -> FloatType {
        value.ln_1p()
    }
    /// Four-quadrant argument of `(x, y)`.
    #[must_use]
    pub fn atan2(y: &FloatType, x: &FloatType) -> FloatType {
        y.atan2(*x)
    }
    /// Approximate real power.
    #[must_use]
    pub fn pow(base: &FloatType, exponent: &FloatType) -> FloatType {
        base.powf(*exponent)
    }
}

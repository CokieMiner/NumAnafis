//! Elementary functions for numeric evaluator values.
//!
//! Arguments are borrowed and results retain the implementing type. Angles
//! are in radians. Scalar real-domain restrictions and ranges are stated below;
//! complex and Clifford implementations extend the domain where representable.
//! Inspection and conversion methods belong to the concrete numeric types.

use core::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, Div, Mul, Neg, Sub},
};

/// Arithmetic and elementary-function contract for evaluator values.
pub trait Numeric:
    Sized
    + Clone
    + Debug
    + Display
    + PartialEq
    + PartialOrd
    + Hash
    + Add<Output = Self>
    + for<'num> Add<&'num Self, Output = Self>
    + Sub<Output = Self>
    + for<'num> Sub<&'num Self, Output = Self>
    + Mul<Output = Self>
    + for<'num> Mul<&'num Self, Output = Self>
    + Div<Output = Self>
    + for<'num> Div<&'num Self, Output = Self>
    + Neg<Output = Self>
{
    // Unary functions

    /// Sine: `sin(self)`.
    #[must_use]
    fn sin(&self) -> Self;

    /// Cosine: `cos(self)`.
    #[must_use]
    fn cos(&self) -> Self;

    /// Tangent: `tan(self) = sin(self)/cos(self)`.
    #[must_use]
    fn tan(&self) -> Self;

    /// Cotangent: `cot(self) = 1/tan(self)`.
    #[must_use]
    fn cot(&self) -> Self;

    /// Secant: `sec(self) = 1/cos(self)`.
    #[must_use]
    fn sec(&self) -> Self;

    /// Cosecant: `csc(self) = 1/sin(self)`.
    #[must_use]
    fn csc(&self) -> Self;

    /// Inverse sine; real inputs in [-1, 1] return values in [-pi/2, pi/2].
    #[must_use]
    fn asin(&self) -> Self;

    /// Inverse cosine; real inputs in [-1, 1] return values in [0, pi].
    #[must_use]
    fn acos(&self) -> Self;

    /// Inverse tangent; real inputs return values in (-pi/2, pi/2).
    #[must_use]
    fn atan(&self) -> Self;

    /// Inverse cotangent: `acot(self)`.
    #[must_use]
    fn acot(&self) -> Self;

    /// Inverse secant: `asec(self)`.
    #[must_use]
    fn asec(&self) -> Self;

    /// Inverse cosecant: `acsc(self)`.
    #[must_use]
    fn acsc(&self) -> Self;

    /// Hyperbolic sine: `sinh(self)`.
    #[must_use]
    fn sinh(&self) -> Self;

    /// Hyperbolic cosine: `cosh(self)`.
    #[must_use]
    fn cosh(&self) -> Self;

    /// Hyperbolic tangent: `tanh(self)`.
    #[must_use]
    fn tanh(&self) -> Self;

    /// Hyperbolic cotangent: `coth(self)`.
    #[must_use]
    fn coth(&self) -> Self;

    /// Hyperbolic secant: `sech(self) = 1/cosh(self)`.
    #[must_use]
    fn sech(&self) -> Self;

    /// Hyperbolic cosecant: `csch(self) = 1/sinh(self)`.
    #[must_use]
    fn csch(&self) -> Self;

    /// Inverse hyperbolic sine: `asinh(self)`.
    #[must_use]
    fn asinh(&self) -> Self;

    /// Inverse hyperbolic cosine; the real-valued domain is [1, infinity).
    #[must_use]
    fn acosh(&self) -> Self;

    /// Inverse hyperbolic tangent; the finite real-valued domain is (-1, 1).
    #[must_use]
    fn atanh(&self) -> Self;

    /// Inverse hyperbolic cotangent; the finite real-valued domain has |x| > 1.
    #[must_use]
    fn acoth(&self) -> Self;

    /// Inverse hyperbolic cosecant: `acsch(self)`.
    #[must_use]
    fn acsch(&self) -> Self;

    /// Inverse hyperbolic secant; the finite real-valued domain is (0, 1].
    #[must_use]
    fn asech(&self) -> Self;

    /// Exponential: `exp(self) = e^self`.
    #[must_use]
    fn exp(&self) -> Self;

    /// Exponential minus one: `expm1(self) = e^self - 1`.
    #[must_use]
    fn expm1(&self) -> Self;

    /// Negative exponential: `exp_neg(self) = e^(-self)`.
    #[must_use]
    fn exp_neg(&self) -> Self;

    /// Natural logarithm; positive real inputs return real values.
    #[must_use]
    fn ln(&self) -> Self;

    /// Natural logarithm of one plus the argument: `log1p(x) = ln(1 + x)`.
    #[must_use]
    fn log1p(&self) -> Self;

    /// Square root: `sqrt(self)`, domain `self >= 0` over the reals.
    #[must_use]
    fn sqrt(&self) -> Self;

    /// Cube root; real inputs retain the real root, including negative inputs.
    #[must_use]
    fn cbrt(&self) -> Self;

    /// Real absolute value or complex scalar modulus.
    #[must_use]
    fn abs(&self) -> Self;

    /// Real sign (-1, 0, or 1), or the direction z/|z| of a nonzero complex scalar.
    #[must_use]
    fn signum(&self) -> Self;

    /// Floor on real values; componentwise floor on complex scalars.
    #[must_use]
    fn floor(&self) -> Self;

    /// Ceiling on real values; componentwise ceiling on complex scalars.
    #[must_use]
    fn ceil(&self) -> Self;

    /// Nearest integer, with halfway cases away from zero; componentwise on complex scalars.
    #[must_use]
    fn round(&self) -> Self;

    /// Fractional part `x - trunc(x)`; componentwise on complex scalars.
    #[must_use]
    fn fract(&self) -> Self;

    /// Arithmetic negation: `-self`.
    #[must_use]
    fn negate(&self) -> Self;

    /// Unnormalized sinc: `sin(self)/self` with `sinc(0) = 1`.
    #[must_use]
    fn sinc(&self) -> Self;
    // Binary / multi-arg functions

    /// Quadrant-aware angle `atan2(y, x)` for real scalars; `self` is y.
    #[must_use]
    fn atan2(&self, x: &Self) -> Self;

    /// Logarithm with explicit base: `log_base(self, base)`.
    #[must_use]
    fn log_base(&self, base: &Self) -> Self;

    /// Power: `self^exp`.
    #[must_use]
    fn pow(&self, exp: &Self) -> Self;

    // TODO: Integrate these methods after the scalar/complex representation and
    // domain, branch, and parameter contracts are verified. Special-function
    // kernels remain in `src/special`; they are not part of Numeric yet.
    // exp_polar also awaits a consistent scalar/Clifford contract.
    // fn erf(&self) -> Self;
    // fn erfc(&self) -> Self;
    // fn gamma(&self) -> Self;
    // fn lgamma(&self) -> Self;
    // fn digamma(&self) -> Self;
    // fn trigamma(&self) -> Self;
    // fn tetragamma(&self) -> Self;
    // fn elliptic_k(&self) -> Self;
    // fn elliptic_e(&self) -> Self;
    // fn zeta(&self) -> Self;
    // fn exp_polar(&self) -> Self;
    // fn besselj(&self, order: &Self) -> Self;
    // fn bessely(&self, order: &Self) -> Self;
    // fn besseli(&self, order: &Self) -> Self;
    // fn besselk(&self, order: &Self) -> Self;
    // fn polygamma(&self, order: &Self) -> Self;
    // fn beta(&self, other: &Self) -> Self;
    // fn zeta_deriv(&self, order: &Self) -> Self;
    // fn lambertw(&self, n: &Self) -> Self;
    // fn hermite(&self, n: &Self) -> Self;
    // fn assoc_legendre(&self, l: &Self, m: &Self) -> Self;
    // fn spherical_harmonic(&self, l: &Self, m: &Self, phi: &Self) -> Self;
}

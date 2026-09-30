//! Number storage for real, complex, and optional Clifford values.

use alloc::boxed::Box;

#[cfg(feature = "clifford")]
use crate::clifford::CliffordNumber;
use crate::error::NumAnafisError;

use super::{Complex, FloatMath, FloatType, IntType, RationalType, Real};

#[derive(Clone, Debug)]
pub enum NumberRepr {
    Real(Real),
    Complex(Box<Complex>),
    #[cfg(feature = "clifford")]
    Clifford(Box<CliffordNumber>),
}

/// One evaluator value; complex components and blade coefficients are real-only.
#[derive(Clone)]
#[non_exhaustive]
pub struct Number(pub(crate) NumberRepr);

impl Number {
    /// Exact native integer.
    #[must_use]
    pub const fn from_int(value: IntType) -> Self {
        Self(NumberRepr::Real(Real::Int(value)))
    }
    /// Reduced rational, normalized to integer when its denominator is one.
    #[must_use]
    pub const fn from_rational(value: RationalType) -> Self {
        Self(NumberRepr::Real(Real::from_rational(value)))
    }
    /// Approximate value; integral floats remain floats.
    #[must_use]
    pub const fn from_float(value: FloatType) -> Self {
        Self(NumberRepr::Real(Real::Float(value)))
    }
    /// Exact zero.
    #[must_use]
    pub const fn zero() -> Self {
        Self::from_int(0)
    }
    /// Exact one.
    #[must_use]
    pub const fn one() -> Self {
        Self::from_int(1)
    }
    /// Exact negative one.
    #[must_use]
    pub const fn neg_one() -> Self {
        Self::from_int(-1)
    }
    /// IEEE NaN.
    #[must_use]
    pub const fn nan() -> Self {
        Self::from_float(FloatMath::nan())
    }
    /// Approximate pi.
    #[must_use]
    pub const fn pi() -> Self {
        Self::from_float(FloatMath::pi())
    }
    /// Approximate value of Euler's number e.
    #[must_use]
    pub const fn e() -> Self {
        Self::from_float(FloatMath::e())
    }
    /// Machine epsilon of the native float backend.
    #[must_use]
    pub const fn epsilon() -> Self {
        Self::from_float(FloatMath::epsilon())
    }
    /// Ordinary scalar imaginary unit.
    #[must_use]
    pub fn i() -> Self {
        Self::from_complex_value(Complex::i())
    }
    /// Constructs `re + i*im`, flattening scalar complex arguments.
    ///
    /// # Panics
    /// Panics when an argument is a general multivector instead of a scalar.
    #[must_use]
    pub fn from_complex(re: Self, im: Self) -> Self {
        if let (NumberRepr::Real(a), NumberRepr::Real(b)) = (&re.0, &im.0) {
            return Self::from_complex_value(Complex::from_parts(a.clone(), b.clone()));
        }

        let a = re
            .into_complex()
            .expect("complex construction requires scalar arguments");
        let b = im
            .into_complex()
            .expect("complex construction requires scalar arguments");
        Self::from_complex_value(Complex::from_parts(a.re - b.im, a.im + b.re))
    }
    pub(crate) fn from_complex_value(value: Complex) -> Self {
        let negative_zero =
            matches!(&value.im, Real::Float(v) if *v == 0.0 && v.is_sign_negative());
        if value.im.is_zero() && !negative_zero {
            if value.im.is_float() {
                value.re.to_float().into()
            } else {
                value.re.into()
            }
        } else {
            Self(NumberRepr::Complex(Box::new(value)))
        }
    }
    pub(crate) fn with_complex(&self, operation: impl FnOnce(&Complex) -> Self) -> Self {
        match &self.0 {
            NumberRepr::Real(value) => operation(&Complex::from_parts(value.clone(), Real::zero())),
            NumberRepr::Complex(value) => operation(value),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(_) => Self::nan(),
        }
    }
    pub(crate) fn into_real(self) -> Result<Real, NumAnafisError> {
        match self.0 {
            NumberRepr::Real(value) => Ok(value),
            NumberRepr::Complex(_) => Err(NumAnafisError::ExpectedReal),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(_) => Err(NumAnafisError::ExpectedReal),
        }
    }
    #[cfg_attr(
        not(feature = "clifford"),
        expect(
            clippy::unnecessary_wraps,
            reason = "Feature-independent checked boundary; Clifford values can fail scalar extraction"
        )
    )]
    pub(crate) fn into_complex(self) -> Result<Complex, NumAnafisError> {
        match self.0 {
            NumberRepr::Real(value) => Ok(Complex::from_parts(value, Real::zero())),
            NumberRepr::Complex(value) => Ok(*value),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(value) => value
                .central_scalar_parts()
                .map(|(re, im)| Complex::from_parts(re.clone(), im.clone()))
                .ok_or(NumAnafisError::ExpectedScalar),
        }
    }
    #[cfg_attr(
        not(feature = "clifford"),
        expect(
            clippy::unnecessary_wraps,
            reason = "Feature-independent inspection contract; general Clifford values have no scalar parts"
        )
    )]
    pub(crate) fn scalar_parts(&self) -> Option<(&Real, &Real)> {
        match &self.0 {
            NumberRepr::Real(value) => Some((value, &Real::Int(0))),
            NumberRepr::Complex(value) => Some((&value.re, &value.im)),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(value) => value.central_scalar_parts(),
        }
    }
    /// Converts a real or complex scalar to typed complex components.
    /// A Clifford value is accepted only when its entire value lies in the
    /// identity/declared-complex-unit subalgebra; its algebra context is discarded.
    ///
    /// # Errors
    /// Rejects values with blades outside the identity and central complex unit.
    pub fn extract_complex(&self) -> Result<Complex, NumAnafisError> {
        self.scalar_parts()
            .map(|(re, im)| Complex::from_parts(re.clone(), im.clone()))
            .ok_or(NumAnafisError::ExpectedScalar)
    }
    /// Borrows the multivector representation.
    #[cfg(feature = "clifford")]
    #[must_use]
    pub fn as_clifford(&self) -> Option<&CliffordNumber> {
        match &self.0 {
            NumberRepr::Clifford(value) => Some(value),
            NumberRepr::Real(_) | NumberRepr::Complex(_) => None,
        }
    }
}

impl From<Real> for Number {
    fn from(component: Real) -> Self {
        Self(NumberRepr::Real(match component {
            Real::Rational(value) => Real::from_rational(value),
            value @ (Real::Int(_) | Real::Float(_)) => value,
        }))
    }
}
impl From<Complex> for Number {
    fn from(value: Complex) -> Self {
        Self::from_complex_value(value)
    }
}
#[cfg(feature = "clifford")]
impl From<CliffordNumber> for Number {
    fn from(value: CliffordNumber) -> Self {
        Self(NumberRepr::Clifford(Box::new(value)))
    }
}

//! Real-only components and direct native backend dispatch.

use alloc::borrow::Cow;

use super::{FloatMath, FloatType, IntType, RationalType};

/// A real component with exact integer and rational tiers.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum Real {
    /// Exact native integer.
    Int(IntType),
    /// Reduced exact rational with unsigned positive denominator.
    Rational(RationalType),
    /// Approximate value, including signed zeros, infinities, and NaNs.
    Float(FloatType),
}

impl From<IntType> for Real {
    fn from(value: IntType) -> Self {
        Self::from_int(value)
    }
}

impl From<RationalType> for Real {
    fn from(value: RationalType) -> Self {
        Self::from_rational(value)
    }
}

impl From<FloatType> for Real {
    fn from(value: FloatType) -> Self {
        Self::from_float(value)
    }
}

// TODO(mp): implement constants, predicates, and explicit primitive conversions
// against owned backend values; native const constructors need adaptation.
impl Real {
    /// IEEE NaN.
    #[must_use]
    pub const fn nan() -> Self {
        Self::Float(FloatMath::nan())
    }
    /// Native machine epsilon.
    #[must_use]
    pub const fn epsilon() -> Self {
        Self::Float(FloatMath::epsilon())
    }
    /// Approximate pi.
    #[must_use]
    pub const fn pi() -> Self {
        Self::Float(FloatMath::pi())
    }
    /// Exact negative one.
    #[must_use]
    pub const fn neg_one() -> Self {
        Self::Int(-1)
    }

    /// Exact zero.
    #[must_use]
    pub const fn zero() -> Self {
        Self::Int(0)
    }
    /// Exact one.
    #[must_use]
    pub const fn one() -> Self {
        Self::Int(1)
    }
    /// Exact native integer.
    #[must_use]
    pub const fn from_int(value: IntType) -> Self {
        Self::Int(value)
    }
    /// Approximate value; integral floats remain floats.
    #[must_use]
    pub const fn from_float(value: FloatType) -> Self {
        Self::Float(value)
    }
    /// Normalizes a denominator-one exact rational.
    #[must_use]
    pub const fn from_rational(value: RationalType) -> Self {
        if value.is_int() {
            Self::Int(*value.numer())
        } else {
            Self::Rational(value)
        }
    }
    /// Whether the approximate representation is active.
    #[must_use]
    pub const fn is_float(&self) -> bool {
        matches!(self, Self::Float(_))
    }
    /// Explicit approximate projection.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        FloatMath::to_f64(&self.float_value())
    }
    /// Converts to the active approximate backend without a primitive projection.
    pub(crate) fn float_value(&self) -> Cow<'_, FloatType> {
        match self {
            Self::Int(value) => Cow::Owned(FloatMath::from_int(value)),
            Self::Rational(value) => {
                Cow::Owned(FloatMath::from_ratio(value.numer(), value.denom()))
            }
            Self::Float(value) => Cow::Borrowed(value),
        }
    }
    /// Explicit conversion to the approximate representation.
    #[must_use]
    pub fn to_float(&self) -> Self {
        Self::Float(self.float_value().into_owned())
    }
    /// Truncating primitive projection with saturating float conversion.
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "Explicit truncating float-to-integer projection"
    )]
    #[must_use]
    pub fn to_int(&self) -> IntType {
        match self {
            Self::Int(value) => *value,
            Self::Rational(value) => value.to_integer(),
            Self::Float(value) => *value as i64,
        }
    }
    /// Whether the represented value is numerically zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match self {
            Self::Int(value) => *value == 0,
            Self::Rational(value) => *value.numer() == 0,
            Self::Float(value) => *value == 0.0,
        }
    }
    /// Whether the represented value is numerically one.
    #[must_use]
    pub fn is_one(&self) -> bool {
        match self {
            Self::Int(value) => *value == 1,
            Self::Rational(value) => value.is_int() && *value.numer() == 1,
            Self::Float(value) => value.total_cmp(&1.0).is_eq(),
        }
    }
    /// Whether the represented value is numerically negative one.
    #[must_use]
    pub fn is_neg_one(&self) -> bool {
        match self {
            Self::Int(value) => *value == -1,
            Self::Rational(value) => value.is_int() && *value.numer() == -1,
            Self::Float(value) => value.total_cmp(&-1.0).is_eq(),
        }
    }
    /// Whether the value is strictly negative; signed zero is excluded.
    #[must_use]
    pub const fn is_negative(&self) -> bool {
        match self {
            Self::Int(value) => *value < 0,
            Self::Rational(value) => *value.numer() < 0,
            Self::Float(value) => *value < 0.0,
        }
    }
    /// Whether the value is strictly positive.
    #[must_use]
    pub fn is_positive(&self) -> bool {
        match self {
            Self::Int(value) => *value > 0,
            Self::Rational(value) => *value.numer() > 0,
            Self::Float(value) => *value > 0.0,
        }
    }
    /// Whether the value is NaN.
    #[must_use]
    pub const fn is_nan(&self) -> bool {
        matches!(self, Self::Float(value) if value.is_nan())
    }
    /// Whether the value is infinite.
    #[must_use]
    pub const fn is_infinite(&self) -> bool {
        matches!(self, Self::Float(value) if value.is_infinite())
    }
    /// Whether the value is finite.
    #[must_use]
    pub const fn is_finite(&self) -> bool {
        !self.is_nan() && !self.is_infinite()
    }
    /// Whether the represented value is integral, including finite integral floats.
    #[must_use]
    pub fn is_int(&self) -> bool {
        match self {
            Self::Int(_) => true,
            Self::Rational(value) => value.is_int(),
            Self::Float(value) => value.is_finite() && FloatMath::fract(value) == 0.0,
        }
    }

    /// Real trunc, preserving the approximate tier for floats.
    #[must_use]
    pub fn trunc(&self) -> Self {
        match self {
            Self::Int(_) => self.clone(),
            Self::Rational(value) => Self::Int(value.trunc().to_integer()),
            Self::Float(value) => Self::Float(FloatMath::trunc(value)),
        }
    }
}

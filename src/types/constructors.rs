//! Convenience constructors for [`Number`].

use core::num::TryFromIntError;

use super::{FloatType, IntMath, Number, RationalMath};

macro_rules! impl_from_exact {
    ($($ty:ty),*) => {
        $(
            impl From<$ty> for Number {
                #[inline]
                fn from(value: $ty) -> Self {
                    Self::from_int(IntMath::from_i64(i64::from(value)))
                }
            }
        )*
    };
}
impl_from_exact!(i8, i16, i32, u8, u16, u32);

impl From<i64> for Number {
    #[inline]
    fn from(value: i64) -> Self {
        Self::from_int(IntMath::from_i64(value))
    }
}

impl TryFrom<u64> for Number {
    type Error = TryFromIntError;

    fn try_from(value: u64) -> Result<Self, Self::Error> {
        i64::try_from(value).map(Self::from_int)
    }
}

impl TryFrom<usize> for Number {
    type Error = TryFromIntError;

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        i64::try_from(value).map(Self::from_int)
    }
}

impl From<f32> for Number {
    #[inline]
    fn from(value: f32) -> Self {
        Self::from_float(FloatType::from(value))
    }
}

impl From<f64> for Number {
    #[inline]
    fn from(value: f64) -> Self {
        Self::from_float(value)
    }
}

/// Create a [`Number`] from any convertible numeric value.
#[inline]
#[must_use]
pub fn n<T: Into<Number>>(val: T) -> Number {
    val.into()
}

/// Create a rational [`Number`] from integer numerator and denominator.
#[inline]
#[must_use]
pub fn r(num: i64, den: i64) -> Number {
    if den == 0 {
        return Number::nan();
    }
    Number::from_rational(RationalMath::new(
        IntMath::from_i64(num),
        IntMath::from_i64(den),
    ))
}

/// Create a complex [`Number`] from real and imaginary parts.
#[inline]
#[must_use]
pub fn c<R: Into<Number>, I: Into<Number>>(re: R, im: I) -> Number {
    Number::from_complex(re.into(), im.into())
}

/// Returns the imaginary unit $i = 0 + 1i$.
#[inline]
#[must_use]
pub fn i() -> Number {
    Number::i()
}

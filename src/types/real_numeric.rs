//! Real-domain evaluator kernels; domain extension belongs to Number.

use crate::traits::Numeric;

use super::{FloatMath, IntMath, RationalType, Real};

impl Numeric for Real {
    fn sin(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::sin(&self.float_value()))
    }
    fn cos(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::one();
        }
        Self::Float(FloatMath::cos(&self.float_value()))
    }
    fn tan(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::tan(&self.float_value()))
    }
    fn sinh(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::sinh(&self.float_value()))
    }
    fn cosh(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::one();
        }
        Self::Float(FloatMath::cosh(&self.float_value()))
    }
    fn tanh(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::tanh(&self.float_value()))
    }
    fn exp(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::one();
        }
        Self::Float(FloatMath::exp(&self.float_value()))
    }
    fn ln(&self) -> Self {
        if !self.is_float() && self.is_one() {
            return Self::zero();
        }
        Self::Float(FloatMath::ln(&self.float_value()))
    }
    fn sqrt(&self) -> Self {
        match self {
            Self::Int(value) => {
                if let Some(root) = IntMath::perfect_square(value) {
                    return Self::Int(root);
                }
            }
            Self::Rational(value) => {
                if let (Some(numer), Some(denom)) = (
                    IntMath::perfect_square(value.numer()),
                    IntMath::perfect_square_unsigned(*value.denom()),
                ) {
                    return Self::from_rational(RationalType::from_parts(numer, denom));
                }
            }
            Self::Float(_) => {}
        }
        Self::Float(FloatMath::sqrt(&self.float_value()))
    }
    fn cbrt(&self) -> Self {
        match self {
            Self::Int(value) => {
                if let Some(root) = IntMath::perfect_cube(value) {
                    return Self::Int(root);
                }
            }
            Self::Rational(value) => {
                if let (Some(numer), Some(denom)) = (
                    IntMath::perfect_cube(value.numer()),
                    IntMath::perfect_cube_unsigned(*value.denom()),
                ) {
                    return Self::from_rational(RationalType::from_parts(numer, denom));
                }
            }
            Self::Float(_) => {}
        }
        Self::Float(FloatMath::cbrt(&self.float_value()))
    }
    fn abs(&self) -> Self {
        match self {
            Self::Int(value) => value.checked_abs().map_or_else(
                || Self::Float(FloatMath::abs(&self.float_value())),
                Self::Int,
            ),
            Self::Rational(value) if self.is_negative() => value.numer().checked_abs().map_or_else(
                || Self::Float(FloatMath::abs(&self.float_value())),
                |numer| Self::from_rational(RationalType::from_parts(numer, *value.denom())),
            ),
            Self::Rational(_) => self.clone(),
            Self::Float(value) => Self::Float(FloatMath::abs(value)),
        }
    }
    fn floor(&self) -> Self {
        match self {
            Self::Int(_) => self.clone(),
            Self::Rational(value) => Self::Int(value.floor().to_integer()),
            Self::Float(value) => Self::Float(FloatMath::floor(value)),
        }
    }
    fn ceil(&self) -> Self {
        match self {
            Self::Int(_) => self.clone(),
            Self::Rational(value) => Self::Int(value.ceil().to_integer()),
            Self::Float(value) => Self::Float(FloatMath::ceil(value)),
        }
    }
    fn round(&self) -> Self {
        match self {
            Self::Int(_) => self.clone(),
            Self::Rational(value) => Self::Int(value.round().to_integer()),
            Self::Float(value) => Self::Float(FloatMath::round(value)),
        }
    }
    fn fract(&self) -> Self {
        match self {
            Self::Int(_) => Self::zero(),
            Self::Rational(value) => Self::from_rational(value.fract()),
            Self::Float(value) => Self::Float(FloatMath::fract(value)),
        }
    }
    fn atan2(&self, x: &Self) -> Self {
        Self::Float(FloatMath::atan2(&self.float_value(), &x.float_value()))
    }
    fn asin(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::asin(&self.float_value()))
    }
    fn acos(&self) -> Self {
        Self::Float(FloatMath::acos(&self.float_value()))
    }
    fn atan(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::atan(&self.float_value()))
    }
    fn asinh(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::asinh(&self.float_value()))
    }
    fn acosh(&self) -> Self {
        Self::Float(FloatMath::acosh(&self.float_value()))
    }
    fn atanh(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::atanh(&self.float_value()))
    }
    fn expm1(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::expm1(&self.float_value()))
    }
    fn log1p(&self) -> Self {
        if !self.is_float() && self.is_zero() {
            return Self::zero();
        }
        Self::Float(FloatMath::log1p(&self.float_value()))
    }
    fn cot(&self) -> Self {
        Self::one() / self.tan()
    }
    fn sec(&self) -> Self {
        Self::one() / self.cos()
    }
    fn csc(&self) -> Self {
        Self::one() / self.sin()
    }
    fn coth(&self) -> Self {
        Self::one() / self.tanh()
    }
    fn sech(&self) -> Self {
        Self::one() / self.cosh()
    }
    fn csch(&self) -> Self {
        Self::one() / self.sinh()
    }
    fn acot(&self) -> Self {
        (Self::one() / self).atan()
    }
    fn asec(&self) -> Self {
        (Self::one() / self).acos()
    }
    fn acsc(&self) -> Self {
        (Self::one() / self).asin()
    }
    fn acoth(&self) -> Self {
        (Self::one() / self).atanh()
    }
    fn acsch(&self) -> Self {
        (Self::one() / self).asinh()
    }
    fn asech(&self) -> Self {
        (Self::one() / self).acosh()
    }
    fn exp_neg(&self) -> Self {
        (-self).exp()
    }
    fn negate(&self) -> Self {
        -self
    }
    fn sinc(&self) -> Self {
        if self.is_zero() {
            Self::one()
        } else {
            self.sin() / self
        }
    }
    fn signum(&self) -> Self {
        if self.is_nan() {
            Self::nan()
        } else if self.is_zero() {
            self.clone()
        } else if self.is_negative() {
            Self::neg_one()
        } else {
            Self::one()
        }
    }
    fn log_base(&self, base: &Self) -> Self {
        {
            let numerator = self.ln();
            numerator / base.ln()
        }
    }
    fn pow(&self, exponent: &Self) -> Self {
        if let Self::Int(e) = exponent
            && !self.is_float()
        {
            return integer_power(self, *e);
        }
        Self::Float(FloatMath::pow(&self.float_value(), &exponent.float_value()))
    }
}

// TODO(mp): traverse exponent magnitude bits directly, without i64 projection.
fn integer_power(base: &Real, exponent: i64) -> Real {
    let mut remaining = exponent.unsigned_abs();
    let mut factor = if exponent < 0 {
        Real::one() / base
    } else {
        base.clone()
    };
    let mut result = Real::one();
    while remaining != 0 {
        if remaining & 1 != 0 {
            result = result * &factor;
        }
        remaining >>= 1;
        if remaining != 0 {
            factor = &factor * &factor;
        }
    }
    result
}

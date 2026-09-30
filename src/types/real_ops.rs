//! Direct real arithmetic with exact rational promotion.

use core::{
    fmt::{Display, Formatter, Result as FmtResult},
    ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign},
};

use super::{FloatMath, IntMath, RationalMath, Real};

impl Add for &Real {
    type Output = Real;
    fn add(self, rhs: Self) -> Real {
        match (self, rhs) {
            (Real::Int(a), Real::Int(b)) => Real::Int(IntMath::add(a, b)),
            (Real::Int(a), Real::Rational(b)) => {
                Real::from_rational(RationalMath::add(&RationalMath::from_integer(*a), b))
            }
            (Real::Rational(a), Real::Int(b)) => {
                Real::from_rational(RationalMath::add(a, &RationalMath::from_integer(*b)))
            }
            (Real::Rational(a), Real::Rational(b)) => Real::from_rational(RationalMath::add(a, b)),
            _ => Real::Float(FloatMath::add(&self.float_value(), &rhs.float_value())),
        }
    }
}
impl Add for Real {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        self += &rhs;
        self
    }
}
impl Add<&Self> for Real {
    type Output = Self;
    fn add(mut self, rhs: &Self) -> Self {
        self += rhs;
        self
    }
}
impl Sub for &Real {
    type Output = Real;
    fn sub(self, rhs: Self) -> Real {
        match (self, rhs) {
            (Real::Int(a), Real::Int(b)) => Real::Int(IntMath::sub(a, b)),
            (Real::Int(a), Real::Rational(b)) => {
                Real::from_rational(RationalMath::sub(&RationalMath::from_integer(*a), b))
            }
            (Real::Rational(a), Real::Int(b)) => {
                Real::from_rational(RationalMath::sub(a, &RationalMath::from_integer(*b)))
            }
            (Real::Rational(a), Real::Rational(b)) => Real::from_rational(RationalMath::sub(a, b)),
            _ => Real::Float(FloatMath::sub(&self.float_value(), &rhs.float_value())),
        }
    }
}
impl Sub for Real {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        self -= &rhs;
        self
    }
}
impl Sub<&Self> for Real {
    type Output = Self;
    fn sub(mut self, rhs: &Self) -> Self {
        self -= rhs;
        self
    }
}
impl Mul for &Real {
    type Output = Real;
    fn mul(self, rhs: Self) -> Real {
        match (self, rhs) {
            (Real::Int(a), Real::Int(b)) => Real::Int(IntMath::mul(a, b)),
            (Real::Int(a), Real::Rational(b)) => {
                Real::from_rational(RationalMath::mul(&RationalMath::from_integer(*a), b))
            }
            (Real::Rational(a), Real::Int(b)) => {
                Real::from_rational(RationalMath::mul(a, &RationalMath::from_integer(*b)))
            }
            (Real::Rational(a), Real::Rational(b)) => Real::from_rational(RationalMath::mul(a, b)),
            _ => Real::Float(FloatMath::mul(&self.float_value(), &rhs.float_value())),
        }
    }
}
impl Mul for Real {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        &self * &rhs
    }
}
impl Mul<&Self> for Real {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self {
        &self * rhs
    }
}
impl Div for &Real {
    type Output = Real;
    fn div(self, rhs: Self) -> Real {
        if rhs.is_zero() {
            return Real::Float(FloatMath::div(&self.float_value(), &rhs.float_value()));
        }
        match (self, rhs) {
            (Real::Int(a), Real::Int(b)) => Real::from_rational(RationalMath::new(*a, *b)),
            (Real::Int(a), Real::Rational(b)) => {
                Real::from_rational(RationalMath::div(&RationalMath::from_integer(*a), b))
            }
            (Real::Rational(a), Real::Int(b)) => {
                Real::from_rational(RationalMath::div(a, &RationalMath::from_integer(*b)))
            }
            (Real::Rational(a), Real::Rational(b)) => Real::from_rational(RationalMath::div(a, b)),
            _ => Real::Float(FloatMath::div(&self.float_value(), &rhs.float_value())),
        }
    }
}
impl Div for Real {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        &self / &rhs
    }
}
impl Div<&Self> for Real {
    type Output = Self;
    fn div(self, rhs: &Self) -> Self {
        &self / rhs
    }
}

impl Neg for &Real {
    type Output = Real;
    fn neg(self) -> Real {
        match self {
            Real::Int(value) => Real::Int(IntMath::wrapping_neg(value)),
            Real::Rational(value) => Real::from_rational(RationalMath::neg(value)),
            Real::Float(value) => Real::Float(FloatMath::neg(value)),
        }
    }
}
impl Neg for Real {
    type Output = Self;
    fn neg(self) -> Self {
        -&self
    }
}
impl AddAssign<&Self> for Real {
    fn add_assign(&mut self, rhs: &Self) {
        match (self, rhs) {
            (Self::Int(lhs), Self::Int(other)) => IntMath::add_assign(lhs, other),
            (Self::Float(lhs), other) => FloatMath::add_assign(lhs, &other.float_value()),
            (target, other) => *target = &*target + other,
        }
    }
}
impl SubAssign<&Self> for Real {
    fn sub_assign(&mut self, rhs: &Self) {
        match (self, rhs) {
            (Self::Int(lhs), Self::Int(other)) => IntMath::sub_assign(lhs, other),
            (Self::Float(lhs), other) => FloatMath::sub_assign(lhs, &other.float_value()),
            (target, other) => *target = &*target - other,
        }
    }
}
impl Display for Real {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            Self::Int(value) => Display::fmt(value, f),
            Self::Rational(value) => Display::fmt(value, f),
            Self::Float(value) => Display::fmt(value, f),
        }
    }
}

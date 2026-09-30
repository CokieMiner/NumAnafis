//! Typed complex elementary kernels using the principal logarithm.

use crate::traits::Numeric;

use super::{Complex, Real};

impl Numeric for Complex {
    fn sin(&self) -> Self {
        Self::sin(self)
    }
    fn cos(&self) -> Self {
        Self::cos(self)
    }
    fn tan(&self) -> Self {
        Self::tan(self)
    }
    fn sinh(&self) -> Self {
        Self::sinh(self)
    }
    fn cosh(&self) -> Self {
        Self::cosh(self)
    }
    fn tanh(&self) -> Self {
        Self::tanh(self)
    }
    fn exp(&self) -> Self {
        Self::exp(self)
    }
    fn ln(&self) -> Self {
        Self::ln(self)
    }
    fn sqrt(&self) -> Self {
        Self::sqrt(self)
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
    fn asin(&self) -> Self {
        -Self::i() * (Self::i() * self + (Self::one() - self * self).sqrt()).ln()
    }
    fn acos(&self) -> Self {
        Self::from_parts(Real::pi() / Real::from_int(2), Real::zero()) - self.asin()
    }
    fn atan(&self) -> Self {
        let iz = Self::i() * self;
        (Self::one() - &iz).ln().sub(&(Self::one() + iz).ln())
            * Self::from_parts(Real::zero(), Real::one() / Real::from_int(2))
    }
    fn asinh(&self) -> Self {
        (self.clone() + (self * self + Self::one()).sqrt()).ln()
    }
    fn acosh(&self) -> Self {
        (self.clone() + (self.clone() + Self::one()).sqrt() * (self.clone() - Self::one()).sqrt())
            .ln()
    }
    fn atanh(&self) -> Self {
        ((Self::one() + self).ln() - (Self::one() - self).ln())
            / Self::from_parts(Real::from_int(2), Real::zero())
    }
    fn expm1(&self) -> Self {
        let half_sin = (&self.im / &Real::from_int(2)).sin();
        Self::from_parts(
            self.re.expm1() * self.im.cos() - Real::from_int(2) * &half_sin * half_sin,
            self.re.exp() * self.im.sin(),
        )
    }
    fn log1p(&self) -> Self {
        if self.im.is_zero() && !(&self.re + &Real::one()).is_negative() {
            return Self::from_parts(self.re.log1p(), self.im.clone());
        }
        if self.abs() <= Real::epsilon().sqrt() {
            let z2 = self * self;
            let z3 = &z2 * self;
            return self.clone() - z2 / Self::from_parts(Real::from_int(2), Real::zero())
                + z3 / Self::from_parts(Real::from_int(3), Real::zero());
        }
        (self.clone() + Self::one()).ln()
    }
    fn exp_neg(&self) -> Self {
        (-self.clone()).exp()
    }
    fn cbrt(&self) -> Self {
        if self.im.is_zero() {
            Self::from_parts(self.re.cbrt(), self.im.clone())
        } else {
            self.pow(&Self::from_parts(
                Real::one() / Real::from_int(3),
                Real::zero(),
            ))
        }
    }
    fn abs(&self) -> Self {
        Self::from_parts(Self::abs(self), Real::zero())
    }
    fn signum(&self) -> Self {
        if self.re.is_zero() && self.im.is_zero() {
            self.clone()
        } else {
            self.clone() / Self::from_parts(Self::abs(self), Real::zero())
        }
    }
    fn floor(&self) -> Self {
        Self::from_parts(self.re.floor(), self.im.floor())
    }
    fn ceil(&self) -> Self {
        Self::from_parts(self.re.ceil(), self.im.ceil())
    }
    fn round(&self) -> Self {
        Self::from_parts(self.re.round(), self.im.round())
    }
    fn fract(&self) -> Self {
        Self::from_parts(self.re.fract(), self.im.fract())
    }
    fn negate(&self) -> Self {
        -self.clone()
    }
    fn sinc(&self) -> Self {
        if self.re.is_zero() && self.im.is_zero() {
            Self::one()
        } else {
            self.sin() / self
        }
    }
    fn atan2(&self, x: &Self) -> Self {
        if self.im.is_zero() && x.im.is_zero() {
            Self::from_parts(self.re.atan2(&x.re), Real::zero())
        } else {
            Self::from_parts(Real::nan(), Real::nan())
        }
    }
    fn log_base(&self, base: &Self) -> Self {
        {
            let numerator = self.ln();
            numerator / base.ln()
        }
    }
    fn pow(&self, exponent: &Self) -> Self {
        if exponent.im.is_zero()
            && let Real::Int(e) = &exponent.re
        {
            let mut remaining = e.unsigned_abs();
            let mut factor = if *e < 0 {
                Self::one() / self
            } else {
                self.clone()
            };
            let mut result = Self::one();
            while remaining != 0 {
                if remaining & 1 != 0 {
                    result = result * &factor;
                }
                remaining >>= 1;
                if remaining != 0 {
                    factor = &factor * &factor;
                }
            }
            return result;
        }
        (exponent * self.ln()).exp()
    }
}

//! Complex number representation with [`Real`] components.

use core::{
    cmp::Ordering,
    fmt::{Debug, Display, Formatter, Result as FmtResult},
    hash::{Hash, Hasher},
    ops::{Add, Div, Mul, Neg, Sub},
};

use crate::{error::NumAnafisError, traits::Numeric};

use super::{FloatMath, IntMath, Number, Real};

/// A complex number $z = a + bi$ composed of two [`Real`] scalars.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct Complex {
    /// The real component $a = \operatorname{Re}(z)$.
    pub re: Real,
    /// The imaginary component $b = \operatorname{Im}(z)$.
    pub im: Real,
}

#[expect(
    clippy::same_name_method,
    reason = "Inherent binary arithmetic methods mirror standard operator traits"
)]
impl Complex {
    /// Creates a new complex number from real and imaginary parts.
    #[inline]
    #[must_use]
    pub(crate) const fn from_parts(re: Real, im: Real) -> Self {
        Self { re, im }
    }

    /// Constructs a typed complex value from two real numeric values.
    /// Signed zeros are preserved, including on branch cuts.
    ///
    /// # Errors
    /// Rejects complex or Clifford components.
    pub fn new(re: impl Into<Number>, im: impl Into<Number>) -> Result<Self, NumAnafisError> {
        Ok(Self::from_parts(
            re.into().into_real()?,
            im.into().into_real()?,
        ))
    }

    /// Creates a complex zero $0 + 0i$.
    #[inline]
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            re: Real::zero(),
            im: Real::zero(),
        }
    }

    /// Creates a complex one $1 + 0i$.
    #[inline]
    #[must_use]
    pub const fn one() -> Self {
        Self {
            re: Real::one(),
            im: Real::zero(),
        }
    }

    /// Creates the imaginary unit $i = 0 + 1i$.
    #[inline]
    #[must_use]
    pub const fn i() -> Self {
        Self {
            re: Real::zero(),
            im: Real::one(),
        }
    }

    /// Returns the complex conjugate $\bar{z} = a - bi$.
    #[inline]
    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: -&self.im,
        }
    }

    /// Returns the squared modulus $|z|^2 = a^2 + b^2$.
    #[inline]
    #[must_use]
    pub fn norm_sq(&self) -> Real {
        &self.re * &self.re + &self.im * &self.im
    }

    /// Returns the modulus $|z| = \sqrt{a^2 + b^2}$.
    #[inline]
    #[must_use]
    pub fn abs(&self) -> Real {
        if self.im.is_zero() {
            return self.re.abs();
        }
        if self.re.is_zero() {
            return self.im.abs();
        }
        if let (Real::Int(re), Real::Int(im)) = (&self.re, &self.im) {
            if let Some(value) = IntMath::perfect_hypot(re, im) {
                return Real::from_int(value);
            }
            return Real::from_float(FloatMath::hypot(
                &self.re.float_value(),
                &self.im.float_value(),
            ));
        }

        if self.re.is_float() || self.im.is_float() {
            Real::from_float(FloatMath::hypot(
                &self.re.float_value(),
                &self.im.float_value(),
            ))
        } else {
            self.norm_sq().sqrt()
        }
    }

    /// Returns the principal argument $\operatorname{Arg}(z) = \operatorname{atan2}(b, a)$.
    #[inline]
    #[must_use]
    pub fn arg(&self) -> Real {
        self.im.atan2(&self.re)
    }

    /// Complex addition: $(a + bi) + (c + di) = (a + c) + (b + d)i$.
    #[inline]
    #[must_use]
    pub fn add(&self, rhs: &Self) -> Self {
        Self {
            re: &self.re + &rhs.re,
            im: &self.im + &rhs.im,
        }
    }

    /// Complex subtraction: $(a + bi) - (c + di) = (a - c) + (b - d)i$.
    #[inline]
    #[must_use]
    pub fn sub(&self, rhs: &Self) -> Self {
        Self {
            re: &self.re - &rhs.re,
            im: &self.im - &rhs.im,
        }
    }

    /// Complex multiplication: $(a + bi)(c + di) = (ac - bd) + (ad + bc)i$.
    #[inline]
    #[must_use]
    pub fn mul(&self, rhs: &Self) -> Self {
        let re = &self.re * &rhs.re - &self.im * &rhs.im;
        let im = &self.re * &rhs.im + &self.im * &rhs.re;
        Self { re, im }
    }

    /// Complex division: $(a + bi) / (c + di) = \frac{(ac + bd) + (bc - ad)i}{c^2 + d^2}$.
    #[inline]
    #[must_use]
    pub fn div(&self, rhs: &Self) -> Self {
        if (self.re.is_float() || self.im.is_float() || rhs.re.is_float() || rhs.im.is_float())
            && let Some((re, im)) = FloatMath::complex_div(
                &self.re.float_value(),
                &self.im.float_value(),
                &rhs.re.float_value(),
                &rhs.im.float_value(),
            )
        {
            return Self::from_parts(Real::from_float(re), Real::from_float(im));
        }
        let d = rhs.norm_sq();
        let re = (&self.re * &rhs.re + &self.im * &rhs.im) / &d;
        let im = (&self.im * &rhs.re - &self.re * &rhs.im) / &d;
        Self { re, im }
    }

    /// Complex exponential $e^z = e^a (\cos b + i \sin b)$.
    #[must_use]
    pub fn exp(&self) -> Self {
        let r = self.re.exp();
        Self {
            re: &r * &self.im.cos(),
            im: &r * &self.im.sin(),
        }
    }

    /// Complex natural logarithm $\ln(z) = \ln|z| + i \operatorname{Arg}(z)$.
    #[must_use]
    pub fn ln(&self) -> Self {
        Self {
            re: self.abs().ln(),
            im: self.arg(),
        }
    }

    /// Complex square root $\sqrt{z} = \sqrt{|z|} \exp(i \operatorname{Arg}(z) / 2)$.
    #[must_use]
    pub fn sqrt(&self) -> Self {
        if self.im.is_zero() {
            let magnitude = self.re.abs();
            let root = magnitude.sqrt();
            return if self.re.is_negative() {
                let signed_root = if FloatMath::is_sign_negative(&self.im.float_value()) {
                    -root
                } else {
                    root
                };
                Self::from_parts(Real::zero(), signed_root)
            } else {
                Self::from_parts(root, self.im.clone())
            };
        }
        let radius = self.abs();
        let two = Real::from_int(2);
        if self.re.is_negative() {
            let magnitude = ((&radius - &self.re) / &two).sqrt();
            let imaginary = if self.im.is_negative() {
                -magnitude
            } else {
                magnitude
            };
            let real = &self.im / &(&two * &imaginary);
            Self::from_parts(real, imaginary)
        } else {
            let real = ((radius + &self.re) / &two).sqrt();
            let imaginary = &self.im / &(&two * &real);
            Self::from_parts(real, imaginary)
        }
    }

    /// Complex sine: $\sin(a + bi) = \sin a \cosh b + i \cos a \sinh b$.
    #[must_use]
    pub fn sin(&self) -> Self {
        Self {
            re: &self.re.sin() * &self.im.cosh(),
            im: &self.re.cos() * &self.im.sinh(),
        }
    }

    /// Complex cosine: $\cos(a + bi) = \cos a \cosh b - i \sin a \sinh b$.
    #[must_use]
    pub fn cos(&self) -> Self {
        Self {
            re: &self.re.cos() * &self.im.cosh(),
            im: -(&self.re.sin() * &self.im.sinh()),
        }
    }

    /// Complex tangent: $\tan(z) = \sin(z) / \cos(z)$.
    #[must_use]
    pub fn tan(&self) -> Self {
        self.sin().div(&self.cos())
    }

    /// Complex hyperbolic sine: $\sinh(z) = \sinh a \cos b + i \cosh a \sin b$.
    #[must_use]
    pub fn sinh(&self) -> Self {
        Self {
            re: &self.re.sinh() * &self.im.cos(),
            im: &self.re.cosh() * &self.im.sin(),
        }
    }

    /// Complex hyperbolic cosine: $\cosh(z) = \cosh a \cos b + i \sinh a \sin b$.
    #[must_use]
    pub fn cosh(&self) -> Self {
        Self {
            re: &self.re.cosh() * &self.im.cos(),
            im: &self.re.sinh() * &self.im.sin(),
        }
    }

    /// Complex hyperbolic tangent: $\tanh(z) = \sinh(z) / \cosh(z)$.
    #[must_use]
    pub fn tanh(&self) -> Self {
        self.sinh().div(&self.cosh())
    }
}

impl PartialEq for Complex {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.re == other.re && self.im == other.im
    }
}

impl PartialOrd for Complex {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.im.is_zero() && other.im.is_zero() {
            self.re.partial_cmp(&other.re)
        } else {
            (self == other).then_some(Ordering::Equal)
        }
    }
}

impl Hash for Complex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.re.hash(state);
        self.im.hash(state);
    }
}

impl Display for Complex {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        if self.re.is_zero() {
            if self.im.is_one() {
                write!(f, "i")
            } else if self.im.is_neg_one() {
                write!(f, "-i")
            } else {
                write!(f, "{}*i", self.im)
            }
        } else if self.im.is_zero() {
            write!(f, "{}", self.re)
        } else if self.im.is_negative() {
            let neg_im = -&self.im;
            if neg_im.is_one() {
                write!(f, "{} - i", self.re)
            } else {
                write!(f, "{} - {}*i", self.re, neg_im)
            }
        } else if self.im.is_one() {
            write!(f, "{} + i", self.re)
        } else {
            write!(f, "{} + {}*i", self.re, self.im)
        }
    }
}

impl Add for Complex {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::add(&self, &rhs)
    }
}

impl Add<&Self> for Complex {
    type Output = Self;
    #[inline]
    fn add(self, rhs: &Self) -> Self {
        Self::add(&self, rhs)
    }
}

impl Sub for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::sub(&self, &rhs)
    }
}

impl Sub<&Self> for Complex {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: &Self) -> Self {
        Self::sub(&self, rhs)
    }
}

impl Mul for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::mul(&self, &rhs)
    }
}

impl Mul<&Self> for Complex {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: &Self) -> Self {
        Self::mul(&self, rhs)
    }
}

impl Div for Complex {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self::div(&self, &rhs)
    }
}

impl Div<&Self> for Complex {
    type Output = Self;
    #[inline]
    fn div(self, rhs: &Self) -> Self {
        Self::div(&self, rhs)
    }
}

impl Neg for Complex {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl Neg for &Complex {
    type Output = Complex;
    #[inline]
    fn neg(self) -> Complex {
        Complex {
            re: -&self.re,
            im: -&self.im,
        }
    }
}

macro_rules! borrowed_arithmetic {
    ($trait:ident,$method:ident) => {
        impl $trait for &Complex {
            type Output = Complex;
            fn $method(self, rhs: Self) -> Complex {
                Complex::$method(self, rhs)
            }
        }
        impl $trait<Complex> for &Complex {
            type Output = Complex;
            fn $method(self, rhs: Complex) -> Complex {
                Complex::$method(self, &rhs)
            }
        }
    };
}
borrowed_arithmetic!(Add, add);
borrowed_arithmetic!(Sub, sub);
borrowed_arithmetic!(Mul, mul);
borrowed_arithmetic!(Div, div);

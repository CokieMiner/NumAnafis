//! Exact native rationals with signed numerators and unsigned denominators.
//!
//! Wide intermediates are reduced before narrowing. Results beyond the chosen
//! backend's capacity fail explicitly; arithmetic never promotes to floats.

use core::{
    cmp::Ordering,
    fmt::{Display, Formatter, Result as FmtResult},
    ops::{Add, Div, Mul, Neg, Sub},
};

use super::IntType;

/// Reduced rational with a signed numerator and strictly positive denominator.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub struct Ratio<N = IntType, D = u64> {
    numer: N,
    denom: D,
}

// TODO(mp): use unsigned backend magnitudes, cross-cancellation, and gcd/divrem
// instead of fixed-width wide intermediates; preserve canonical borrowed fields.
/// Exact rational representation of the native backend.
pub type RationalType = Ratio<IntType, u64>;

/// Backend facade for rational construction and arithmetic.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct RationalMath;

impl Ratio<IntType, u64> {
    /// Constructs a fraction, accepting either denominator sign.
    ///
    /// # Panics
    /// Panics for denominator zero or an unrepresentable reduced numerator.
    #[must_use]
    pub fn new(numer: IntType, denom: IntType) -> Self {
        Self::from_magnitudes(
            numer.is_negative() != denom.is_negative(),
            u128::from(numer.unsigned_abs()),
            u128::from(denom.unsigned_abs()),
        )
    }
    /// Constructs a fraction with an unsigned denominator.
    ///
    /// # Panics
    /// Panics for denominator zero.
    #[must_use]
    pub fn from_parts(numer: IntType, denom: u64) -> Self {
        Self::from_magnitudes(
            numer.is_negative(),
            u128::from(numer.unsigned_abs()),
            u128::from(denom),
        )
    }
    /// Constructs an exact integer fraction.
    #[must_use]
    pub const fn from_integer(numer: IntType) -> Self {
        Self { numer, denom: 1 }
    }
    /// Signed numerator.
    #[must_use]
    pub const fn numer(&self) -> &IntType {
        &self.numer
    }
    /// Positive unsigned denominator.
    #[must_use]
    pub const fn denom(&self) -> &u64 {
        &self.denom
    }
    /// Whether the reduced fraction represents an integer.
    #[must_use]
    pub const fn is_int(&self) -> bool {
        self.denom == 1
    }
    /// Truncates the quotient toward zero.
    #[expect(
        clippy::integer_division,
        reason = "This explicit projection truncates a rational toward zero"
    )]
    #[must_use]
    pub fn to_integer(&self) -> IntType {
        i64::try_from(i128::from(self.numer) / i128::from(self.denom))
            .expect("quotient cannot grow")
    }
    /// Truncates to an exact integer fraction.
    #[must_use]
    pub fn trunc(&self) -> Self {
        Self::from_integer(self.to_integer())
    }
    /// Exact signed fractional part.
    #[must_use]
    pub fn fract(&self) -> Self {
        let rem = i128::from(self.numer) % i128::from(self.denom);
        Self::from_parts(
            i64::try_from(rem).expect("remainder magnitude is bounded by numerator"),
            self.denom,
        )
    }
    /// Floor as an exact integer fraction.
    #[must_use]
    pub fn floor(&self) -> Self {
        let q = self.to_integer();
        Self::from_integer(if self.numer < 0 && !self.fract().numer.eq(&0) {
            q - 1
        } else {
            q
        })
    }
    /// Ceiling as an exact integer fraction.
    #[must_use]
    pub fn ceil(&self) -> Self {
        let q = self.to_integer();
        Self::from_integer(if self.numer > 0 && !self.fract().numer.eq(&0) {
            q + 1
        } else {
            q
        })
    }
    /// Nearest integer, ties away from zero, without float projection.
    #[must_use]
    pub fn round(&self) -> Self {
        let q = self.to_integer();
        let rem = self.fract().numer.unsigned_abs();
        let increment = u128::from(rem) * 2 >= u128::from(self.denom);
        Self::from_integer(if increment {
            q + self.numer.signum()
        } else {
            q
        })
    }
    #[expect(
        clippy::integer_division,
        reason = "The gcd divides both magnitudes exactly during canonical reduction"
    )]
    fn from_magnitudes(negative: bool, numer: u128, denom: u128) -> Self {
        assert!(denom != 0, "rational denominator must be nonzero");
        let common = gcd(numer, denom);
        let magnitude = numer / common;
        let reduced_denom =
            u64::try_from(denom / common).expect("reduced denominator exceeds native capacity");
        let reduced_numer = if negative && magnitude == (1_u128 << 63) {
            i64::MIN
        } else {
            let signed =
                i64::try_from(magnitude).expect("reduced numerator exceeds native capacity");
            if negative { -signed } else { signed }
        };
        Self {
            numer: reduced_numer,
            denom: reduced_denom,
        }
    }
}

impl Add for &RationalType {
    type Output = RationalType;
    fn add(self, rhs: Self) -> RationalType {
        add_signed(self, rhs, false)
    }
}
impl Sub for &RationalType {
    type Output = RationalType;
    fn sub(self, rhs: Self) -> RationalType {
        add_signed(self, rhs, true)
    }
}
impl Mul for &RationalType {
    type Output = RationalType;
    fn mul(self, rhs: Self) -> RationalType {
        Ratio::from_magnitudes(
            self.numer.is_negative() != rhs.numer.is_negative(),
            u128::from(self.numer.unsigned_abs()) * u128::from(rhs.numer.unsigned_abs()),
            u128::from(self.denom) * u128::from(rhs.denom),
        )
    }
}
impl Div for &RationalType {
    type Output = RationalType;
    fn div(self, rhs: Self) -> RationalType {
        Ratio::from_magnitudes(
            self.numer.is_negative() != rhs.numer.is_negative(),
            u128::from(self.numer.unsigned_abs()) * u128::from(rhs.denom),
            u128::from(self.denom) * u128::from(rhs.numer.unsigned_abs()),
        )
    }
}
impl Neg for &RationalType {
    type Output = RationalType;
    fn neg(self) -> RationalType {
        Ratio::from_magnitudes(
            !self.numer.is_negative(),
            u128::from(self.numer.unsigned_abs()),
            u128::from(self.denom),
        )
    }
}
impl Add for RationalType {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        &self + &rhs
    }
}
impl Sub for RationalType {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        &self - &rhs
    }
}
impl Mul for RationalType {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        &self * &rhs
    }
}
impl Div for RationalType {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        &self / &rhs
    }
}

impl Neg for RationalType {
    type Output = Self;
    fn neg(self) -> Self {
        -&self
    }
}
impl PartialOrd for RationalType {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for RationalType {
    fn cmp(&self, other: &Self) -> Ordering {
        let a = i128::from(self.numer) * i128::from(other.denom);
        let b = i128::from(other.numer) * i128::from(self.denom);
        a.cmp(&b)
    }
}
impl Display for RationalType {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        if self.denom == 1 {
            write!(f, "{}", self.numer)
        } else {
            write!(f, "{}/{}", self.numer, self.denom)
        }
    }
}

impl RationalMath {
    /// Exact integer fraction.
    #[must_use]
    pub const fn from_integer(value: IntType) -> RationalType {
        Ratio::from_integer(value)
    }
    /// Constructs a reduced fraction with a signed denominator.
    ///
    /// # Panics
    /// Panics for denominator zero or a reduced numerator beyond native capacity.
    #[expect(
        clippy::new_ret_no_self,
        reason = "RationalMath is the backend namespace"
    )]
    #[must_use]
    pub fn new(numer: IntType, denom: IntType) -> RationalType {
        Ratio::new(numer, denom)
    }

    /// Exact rational add; wide intermediates are reduced before narrowing.
    ///
    /// # Panics
    /// Panics if the reduced result exceeds native capacity or divides by zero.
    #[must_use]
    pub fn add(lhs: &RationalType, rhs: &RationalType) -> RationalType {
        lhs + rhs
    }
    /// Exact rational sub; wide intermediates are reduced before narrowing.
    ///
    /// # Panics
    /// Panics if the reduced result exceeds native capacity or divides by zero.
    #[must_use]
    pub fn sub(lhs: &RationalType, rhs: &RationalType) -> RationalType {
        lhs - rhs
    }
    /// Exact rational mul; wide intermediates are reduced before narrowing.
    ///
    /// # Panics
    /// Panics if the reduced result exceeds native capacity or divides by zero.
    #[must_use]
    pub fn mul(lhs: &RationalType, rhs: &RationalType) -> RationalType {
        lhs * rhs
    }
    /// Exact rational div; wide intermediates are reduced before narrowing.
    ///
    /// # Panics
    /// Panics if the reduced result exceeds native capacity or divides by zero.
    #[must_use]
    pub fn div(lhs: &RationalType, rhs: &RationalType) -> RationalType {
        lhs / rhs
    }
    /// Exact rational negation.
    ///
    /// # Panics
    /// Panics if the positive numerator is unrepresentable.
    #[must_use]
    pub fn neg(value: &RationalType) -> RationalType {
        -value
    }
    /// Exact comparison without intermediate overflow.
    #[must_use]
    pub fn cmp(lhs: &RationalType, rhs: &RationalType) -> Ordering {
        lhs.cmp(rhs)
    }
}

// Each cross-product fits u128; their sum is at most 2^128 - 2^64.
fn add_signed(lhs: &RationalType, rhs: &RationalType, subtract: bool) -> RationalType {
    let a = u128::from(lhs.numer.unsigned_abs()) * u128::from(rhs.denom);
    let b = u128::from(rhs.numer.unsigned_abs()) * u128::from(lhs.denom);
    let a_negative = lhs.numer.is_negative();
    let b_negative = rhs.numer.is_negative() != subtract;
    let (negative, magnitude) = if a_negative == b_negative {
        (a_negative, a + b)
    } else if a >= b {
        (a_negative, a - b)
    } else {
        (b_negative, b - a)
    };
    Ratio::from_magnitudes(
        negative,
        magnitude,
        u128::from(lhs.denom) * u128::from(rhs.denom),
    )
}

const fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
}

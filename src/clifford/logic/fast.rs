use core::{
    array,
    cmp::Ordering,
    fmt::{Display, Formatter, Result as FmtResult},
    hash::{Hash, Hasher},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use alloc::vec::Vec;

use crate::{error::NumAnafisError, traits::Numeric, types::Real};

use super::{CliffordNumber, GeneratorSet, INLINE_COEFF_COUNT};

/// A fixed-signature Clifford multivector with inline real coefficients.
/// Geometric products use static signature tables. Division and elementary
/// functions delegate to the dynamic implementation.
///
/// `P`, `Q`, `R` define the signature (positive, negative, zero metrics).
/// Requires `P + Q + R <= 5`.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct FastClifford<const P: usize, const Q: usize, const R: usize> {
    /// The array of multivector coefficients, strictly bounded to `INLINE_COEFF_COUNT`.
    pub(crate) coeffs: [Real; INLINE_COEFF_COUNT],
}

impl<const P: usize, const Q: usize, const R: usize> FastClifford<P, Q, R> {
    const CHECK_SIZE: () = assert!(
        P + Q + R <= 5,
        "FastClifford requires P + Q + R <= 5 for inline storage"
    );

    // Signature-only signs are computed at compile time, avoiding dynamic
    // GeneratorSet/Cayley allocation for fixed-size geometric multiplication.
    const PRODUCT_SIGNS: [[i8; 32]; 32] = {
        let mut signs = [[0_i8; 32]; 32];
        let count = 1 << (P + Q + R);
        let mut a = 0_usize;
        while a < count {
            let mut b = 0;
            while b < count {
                let mut factor = 1_i8;
                let mut generator = 0;
                while generator < P + Q + R {
                    if b & (1 << generator) != 0 && ((a >> (generator + 1)).count_ones() & 1) != 0 {
                        factor = -factor;
                    }
                    if a & b & (1 << generator) != 0 {
                        if generator >= P + Q {
                            factor = 0;
                        } else if generator >= P {
                            factor = -factor;
                        }
                    }
                    generator += 1;
                }
                signs[a][b] = factor;
                b += 1;
            }
            a += 1;
        }
        signs
    };
    /// Borrows the active blade coefficients; padding is never exposed.
    #[must_use]
    pub fn coeffs_slice(&self) -> &[Real] {
        let () = Self::CHECK_SIZE;
        &self.coeffs[..1 << Self::active_generators()]
    }
    /// Mutates active real coefficients while preserving zero padding.
    pub fn coeffs_mut_slice(&mut self) -> &mut [Real] {
        let () = Self::CHECK_SIZE;
        &mut self.coeffs[..1 << Self::active_generators()]
    }

    /// Creates a multivector initialized completely to zero.
    #[must_use]
    pub fn zero() -> Self {
        let () = Self::CHECK_SIZE;
        Self {
            coeffs: array::from_fn(|_| Real::zero()),
        }
    }
    /// Returns the total number of active generators mapped to this algebra (`P + Q + R`).
    #[must_use]
    pub const fn active_generators() -> usize {
        P + Q + R
    }
    /// Compares with total ordering over active blade coefficients.
    ///
    /// This coefficient order is a sorting convention.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        let len = 1 << Self::active_generators();
        for (i, c) in self.coeffs.iter().enumerate().take(len) {
            let cmp = c.total_cmp(&other.coeffs[i]);
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }
}

impl<const P: usize, const Q: usize, const R: usize> From<&FastClifford<P, Q, R>>
    for CliffordNumber
{
    fn from(val: &FastClifford<P, Q, R>) -> Self {
        let () = FastClifford::<P, Q, R>::CHECK_SIZE;
        let n = P + Q + R;
        let mut entries = Vec::with_capacity(n);
        for i in 0..P {
            entries.push((u32::try_from(i).expect("i < P+Q+R <= 5 fits in u32"), 1));
        }
        for i in P..P + Q {
            entries.push((u32::try_from(i).expect("i < P+Q+R <= 5 fits in u32"), -1));
        }
        for i in P + Q..n {
            entries.push((u32::try_from(i).expect("i < P+Q+R <= 5 fits in u32"), 0));
        }
        let gens = GeneratorSet::from_validated(entries);

        let mut mv = Self::zero_unchecked(gens);
        let len = 1 << n;
        for (i, coeff) in val.coeffs.iter().enumerate().take(len) {
            mv.coeffs_mut_slice()[i] = coeff.clone();
        }
        mv
    }
}

impl<const P: usize, const Q: usize, const R: usize> FastClifford<P, Q, R> {
    /// Creates a `FastClifford` from a `CliffordNumber` without checking generator match.
    #[must_use]
    fn from_unchecked(val: &CliffordNumber) -> Self {
        let () = Self::CHECK_SIZE;
        let mut coeffs = array::from_fn(|_| Real::zero());
        let len = 1 << (P + Q + R);
        for (i, coeff) in coeffs.iter_mut().enumerate().take(len) {
            *coeff = val.coeff(i).clone();
        }
        Self { coeffs }
    }
}

impl<const P: usize, const Q: usize, const R: usize> TryFrom<&CliffordNumber>
    for FastClifford<P, Q, R>
{
    type Error = NumAnafisError;

    fn try_from(val: &CliffordNumber) -> Result<Self, Self::Error> {
        let () = Self::CHECK_SIZE;
        let gens = val.generator_set();
        if gens.len() != P + Q + R {
            return Err(NumAnafisError::MismatchedGeneratorSet);
        }
        for i in 0..P {
            if gens.metric_at(i) != 1 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        for i in P..P + Q {
            if gens.metric_at(i) != -1 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        for i in P + Q..P + Q + R {
            if gens.metric_at(i) != 0 {
                return Err(NumAnafisError::MismatchedGeneratorSet);
            }
        }
        Ok(Self::from_unchecked(val))
    }
}

impl<const P: usize, const Q: usize, const R: usize> Display for FastClifford<P, Q, R> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        let mv = CliffordNumber::from(self);
        write!(f, "{mv}")
    }
}

impl<const P: usize, const Q: usize, const R: usize> PartialEq for FastClifford<P, Q, R> {
    fn eq(&self, other: &Self) -> bool {
        let len = 1 << Self::active_generators();
        self.coeffs
            .iter()
            .take(len)
            .zip(other.coeffs.iter().take(len))
            .all(|(a, b)| a == b)
    }
}

impl<const P: usize, const Q: usize, const R: usize> PartialOrd for FastClifford<P, Q, R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let len = 1 << Self::active_generators();
        if self.coeffs[1..len].iter().any(|value| !value.is_zero())
            || other.coeffs[1..len].iter().any(|value| !value.is_zero())
        {
            return (self == other).then_some(Ordering::Equal);
        }
        self.coeffs[0].partial_cmp(&other.coeffs[0])
    }
}

impl<const P: usize, const Q: usize, const R: usize> Hash for FastClifford<P, Q, R> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        P.hash(state);
        Q.hash(state);
        R.hash(state);
        let len = 1 << Self::active_generators();
        for c in self.coeffs.iter().take(len) {
            c.hash(state);
        }
    }
}

// Direct zero-allocation implementations for simple coefficient-wise operations.
impl<const P: usize, const Q: usize, const R: usize> Add for FastClifford<P, Q, R> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] += &rhs.coeffs[i];
        }
        self
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Add<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn add(mut self, rhs: &'num Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] += &rhs.coeffs[i];
        }
        self
    }
}

impl<const P: usize, const Q: usize, const R: usize> Sub for FastClifford<P, Q, R> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] -= &rhs.coeffs[i];
        }
        self
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Sub<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn sub(mut self, rhs: &'num Self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] -= &rhs.coeffs[i];
        }
        self
    }
}

// Geometric multiplication uses signature constants; division uses inversion.
impl<const P: usize, const Q: usize, const R: usize> Mul for FastClifford<P, Q, R> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self.product(&rhs)
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Mul<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn mul(self, rhs: &'num Self) -> Self::Output {
        self.product(rhs)
    }
}

impl<const P: usize, const Q: usize, const R: usize> Div for FastClifford<P, Q, R> {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let mv1 = CliffordNumber::from(&self);
        let mv2 = CliffordNumber::from(&rhs);
        let res = mv1.div(&mv2);
        Self::from_unchecked(&res)
    }
}

impl<'num, const P: usize, const Q: usize, const R: usize> Div<&'num Self>
    for FastClifford<P, Q, R>
{
    type Output = Self;
    fn div(self, rhs: &'num Self) -> Self::Output {
        let mv1 = CliffordNumber::from(&self);
        let mv2 = CliffordNumber::from(rhs);
        let res = mv1.div(&mv2);
        Self::from_unchecked(&res)
    }
}

impl<const P: usize, const Q: usize, const R: usize> Neg for FastClifford<P, Q, R> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        let len = 1 << Self::active_generators();
        for coeff in self.coeffs.iter_mut().take(len) {
            *coeff = -&*coeff;
        }
        self
    }
}
// Compound assignment operators

impl<const P: usize, const Q: usize, const R: usize> AddAssign for FastClifford<P, Q, R> {
    fn add_assign(&mut self, rhs: Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] += &rhs.coeffs[i];
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> AddAssign<&Self> for FastClifford<P, Q, R> {
    fn add_assign(&mut self, rhs: &Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] += &rhs.coeffs[i];
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> SubAssign for FastClifford<P, Q, R> {
    fn sub_assign(&mut self, rhs: Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] -= &rhs.coeffs[i];
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> SubAssign<&Self> for FastClifford<P, Q, R> {
    fn sub_assign(&mut self, rhs: &Self) {
        let len = 1 << Self::active_generators();
        for i in 0..len {
            self.coeffs[i] -= &rhs.coeffs[i];
        }
    }
}

impl<const P: usize, const Q: usize, const R: usize> MulAssign for FastClifford<P, Q, R> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.product(&rhs);
    }
}

impl<const P: usize, const Q: usize, const R: usize> MulAssign<&Self> for FastClifford<P, Q, R> {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = self.product(rhs);
    }
}

impl<const P: usize, const Q: usize, const R: usize> DivAssign for FastClifford<P, Q, R> {
    fn div_assign(&mut self, rhs: Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(&rhs);
        *self = Self::from_unchecked(&(mv1 / &mv2));
    }
}

impl<const P: usize, const Q: usize, const R: usize> DivAssign<&Self> for FastClifford<P, Q, R> {
    fn div_assign(&mut self, rhs: &Self) {
        let mv1 = CliffordNumber::from(&*self);
        let mv2 = CliffordNumber::from(rhs);
        *self = Self::from_unchecked(&(mv1 / &mv2));
    }
}

/// Evaluates unary Numeric functions through the dynamic multivector implementation.
macro_rules! delegate_number_unary {
    ($($method:ident),*) => {
        $(
            fn $method(&self) -> Self {
                let mv = CliffordNumber::from(self);
                let res = mv.$method();
                Self::from_unchecked(&res)
            }
        )*
    };
}

macro_rules! delegate_number_binary {
    ($($method:ident),*) => {
        $(
            fn $method(&self, other: &Self) -> Self {
                let mv1 = CliffordNumber::from(self);
                let mv2 = CliffordNumber::from(other);
                let res = mv1.$method(&mv2);
                Self::from_unchecked(&res)
            }
        )*
    };
}

impl<const P: usize, const Q: usize, const R: usize> Numeric for FastClifford<P, Q, R> {
    delegate_number_unary!(
        sin, cos, tan, cot, sec, csc, asin, acos, atan, acot, asec, acsc, sinh, cosh, tanh, coth,
        sech, csch, asinh, acosh, atanh, acoth, acsch, asech, exp, expm1, exp_neg, ln, log1p, sqrt,
        cbrt, abs, signum, floor, ceil, round, fract, negate, sinc
    );

    delegate_number_binary!(atan2, log_base, pow);
}

impl<const P: usize, const Q: usize, const R: usize> FastClifford<P, Q, R> {
    fn product(&self, rhs: &Self) -> Self {
        let mut result = Self::zero();
        let count = 1 << Self::active_generators();
        for a in 0..count {
            if self.coeffs[a].is_zero() {
                continue;
            }
            for b in 0..count {
                let factor = Self::PRODUCT_SIGNS[a][b];
                if factor == 0 || rhs.coeffs[b].is_zero() {
                    continue;
                }
                let raw = &self.coeffs[a] * &rhs.coeffs[b];
                let term = if factor < 0 { -raw } else { raw };
                let target = &mut result.coeffs[a ^ b];
                *target += &term;
            }
        }
        result
    }
}

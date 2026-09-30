use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::types::{FloatType, IntType, Real};

use super::CliffordNumber;
// Neg

impl Neg for CliffordNumber {
    type Output = Self;
    fn neg(self) -> Self {
        let mut out = self;
        for value in out.coeffs_mut_slice() {
            *value = -&*value;
        }
        out
    }
}

impl Neg for &CliffordNumber {
    type Output = CliffordNumber;
    fn neg(self) -> CliffordNumber {
        -(self.clone())
    }
}
// Add

impl Add for CliffordNumber {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut out = self;
        out += &rhs;
        out
    }
}

impl Add<&Self> for CliffordNumber {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self {
        let mut out = self;
        out += rhs;
        out
    }
}

impl Add<CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn add(self, rhs: CliffordNumber) -> CliffordNumber {
        self.clone() + rhs
    }
}

impl Add<&CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn add(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.clone() + rhs
    }
}
// Sub

impl Sub for CliffordNumber {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut out = self;
        out -= &rhs;
        out
    }
}

impl Sub<&Self> for CliffordNumber {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self {
        let mut out = self;
        out -= rhs;
        out
    }
}

impl Sub<CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn sub(self, rhs: CliffordNumber) -> CliffordNumber {
        self.clone() - rhs
    }
}

impl Sub<&CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn sub(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.clone() - rhs
    }
}
// Mul (geometric product)

impl Mul for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.geometric_mul(&rhs)
    }
}

impl Mul<&Self> for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self {
        self.geometric_mul(rhs)
    }
}

impl Mul<CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: CliffordNumber) -> CliffordNumber {
        self.geometric_mul(&rhs)
    }
}

impl Mul<&CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.geometric_mul(rhs)
    }
}
// Div (right-division: a / b = a·b⁻¹)
//
// In a non-commutative algebra `a / b` is ambiguous — this implements
// **right-division** (`a * geometric_inverse(b)`), matching the convention
// `A/B = AB⁻¹` from linear algebra.  Use `geometric_inverse` + `geometric_mul`
// directly for left-division (`A⁻¹B`).

impl Div for CliffordNumber {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self.geometric_mul(&rhs.geometric_inverse())
    }
}

impl Div<&Self> for CliffordNumber {
    type Output = Self;
    fn div(self, rhs: &Self) -> Self {
        self.geometric_mul(&rhs.geometric_inverse())
    }
}

impl Div<CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn div(self, rhs: CliffordNumber) -> CliffordNumber {
        self.geometric_mul(&rhs.geometric_inverse())
    }
}

impl Div<&CliffordNumber> for &CliffordNumber {
    type Output = CliffordNumber;
    fn div(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.geometric_mul(&rhs.geometric_inverse())
    }
}
// Real interop — CliffordNumber × Real

impl Mul<Real> for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: Real) -> Self {
        scale(&self, &rhs)
    }
}

impl Mul<&Real> for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: &Real) -> Self {
        scale(&self, rhs)
    }
}

impl Mul<Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: Real) -> CliffordNumber {
        scale(self, &rhs)
    }
}

impl Mul<&Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: &Real) -> CliffordNumber {
        scale(self, rhs)
    }
}

// Real × CliffordNumber
impl Mul<CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn mul(self, rhs: CliffordNumber) -> CliffordNumber {
        scale(&rhs, &self)
    }
}

impl Mul<&CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn mul(self, rhs: &CliffordNumber) -> CliffordNumber {
        scale(rhs, &self)
    }
}

impl Mul<CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn mul(self, rhs: CliffordNumber) -> CliffordNumber {
        scale(&rhs, self)
    }
}

impl Mul<&CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn mul(self, rhs: &CliffordNumber) -> CliffordNumber {
        scale(rhs, self)
    }
}

// CliffordNumber × IntType
impl Mul<IntType> for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: IntType) -> Self {
        scale(&self, &Real::Int(rhs))
    }
}

impl Mul<IntType> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: IntType) -> CliffordNumber {
        scale(self, &Real::Int(rhs))
    }
}

// CliffordNumber × FloatType
impl Mul<FloatType> for CliffordNumber {
    type Output = Self;
    fn mul(self, rhs: FloatType) -> Self {
        scale(&self, &Real::from_float(rhs))
    }
}

impl Mul<FloatType> for &CliffordNumber {
    type Output = CliffordNumber;
    fn mul(self, rhs: FloatType) -> CliffordNumber {
        scale(self, &Real::from_float(rhs))
    }
}
// Real + CliffordNumber

impl Add<CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn add(self, rhs: CliffordNumber) -> CliffordNumber {
        scalar_plus_mv(&self, rhs)
    }
}

impl Add<&CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn add(self, rhs: &CliffordNumber) -> CliffordNumber {
        self + rhs.clone()
    }
}

impl Add<CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn add(self, rhs: CliffordNumber) -> CliffordNumber {
        self.clone() + rhs
    }
}

impl Add<&CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn add(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.clone() + rhs.clone()
    }
}
// CliffordNumber + Real

impl Add<Real> for CliffordNumber {
    type Output = Self;
    fn add(self, rhs: Real) -> Self {
        scalar_plus_mv(&rhs, self)
    }
}

impl Add<&Real> for CliffordNumber {
    type Output = Self;
    fn add(self, rhs: &Real) -> Self {
        scalar_plus_mv(rhs, self)
    }
}

impl Add<Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn add(self, rhs: Real) -> CliffordNumber {
        self.clone() + rhs
    }
}

impl Add<&Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn add(self, rhs: &Real) -> CliffordNumber {
        self.clone() + rhs.clone()
    }
}
// Real - CliffordNumber

impl Sub<CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn sub(self, rhs: CliffordNumber) -> CliffordNumber {
        self + (-rhs)
    }
}

impl Sub<&CliffordNumber> for Real {
    type Output = CliffordNumber;
    fn sub(self, rhs: &CliffordNumber) -> CliffordNumber {
        self - rhs.clone()
    }
}

impl Sub<CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn sub(self, rhs: CliffordNumber) -> CliffordNumber {
        self.clone() - rhs
    }
}

impl Sub<&CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn sub(self, rhs: &CliffordNumber) -> CliffordNumber {
        self.clone() - rhs.clone()
    }
}
// CliffordNumber - Real

impl Sub<Real> for CliffordNumber {
    type Output = Self;
    fn sub(self, rhs: Real) -> Self {
        self + (-rhs)
    }
}

impl Sub<&Real> for CliffordNumber {
    type Output = Self;
    fn sub(self, rhs: &Real) -> Self {
        self - rhs.clone()
    }
}

impl Sub<Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn sub(self, rhs: Real) -> CliffordNumber {
        self.clone() - rhs
    }
}

impl Sub<&Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn sub(self, rhs: &Real) -> CliffordNumber {
        self.clone() - rhs.clone()
    }
}
// Compound assignment operators

impl AddAssign for CliffordNumber {
    fn add_assign(&mut self, rhs: Self) {
        if self.generator_set() == rhs.generator_set() {
            let limit = self.blade_count();
            for i in 0..limit {
                self.coeffs_mut_slice()[i] += &rhs.coeffs_slice()[i];
            }
        } else {
            *self = self.clone() + rhs;
        }
    }
}

impl AddAssign<&Self> for CliffordNumber {
    fn add_assign(&mut self, rhs: &Self) {
        if self.generator_set() == rhs.generator_set() {
            let limit = self.blade_count();
            for i in 0..limit {
                self.coeffs_mut_slice()[i] += &rhs.coeffs_slice()[i];
            }
        } else {
            *self = self.clone() + rhs.clone();
        }
    }
}

impl SubAssign for CliffordNumber {
    fn sub_assign(&mut self, rhs: Self) {
        if self.generator_set() == rhs.generator_set() {
            let limit = self.blade_count();
            for i in 0..limit {
                self.coeffs_mut_slice()[i] -= &rhs.coeffs_slice()[i];
            }
        } else {
            *self = self.clone() - rhs;
        }
    }
}

impl SubAssign<&Self> for CliffordNumber {
    fn sub_assign(&mut self, rhs: &Self) {
        if self.generator_set() == rhs.generator_set() {
            let limit = self.blade_count();
            for i in 0..limit {
                self.coeffs_mut_slice()[i] -= &rhs.coeffs_slice()[i];
            }
        } else {
            *self = self.clone() - rhs.clone();
        }
    }
}

impl MulAssign for CliffordNumber {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.geometric_mul(&rhs);
    }
}

impl MulAssign<&Self> for CliffordNumber {
    fn mul_assign(&mut self, rhs: &Self) {
        *self = self.geometric_mul(rhs);
    }
}

impl DivAssign for CliffordNumber {
    fn div_assign(&mut self, rhs: Self) {
        *self = self.geometric_mul(&rhs.geometric_inverse());
    }
}

impl DivAssign<&Self> for CliffordNumber {
    fn div_assign(&mut self, rhs: &Self) {
        *self = self.geometric_mul(&rhs.geometric_inverse());
    }
}
// Shared leaf kernels
//
// A grade-0 scalar commutes with every blade, so scalar-multivector arithmetic
// never needs the geometric product: scaling touches each blade coefficient
// once, and scalar addition touches only the scalar blade (index 0).

fn scale(mv: &CliffordNumber, scalar: &Real) -> CliffordNumber {
    let mut out = mv.clone();
    let limit = mv.blade_count();
    for i in 0..limit {
        out.coeffs_mut_slice()[i] = &mv.coeffs_slice()[i] * scalar;
    }
    out
}

fn scalar_plus_mv(scalar: &Real, mut mv: CliffordNumber) -> CliffordNumber {
    let new_scalar = &mv.coeffs_slice()[0] + scalar;
    mv.coeffs_mut_slice()[0] = new_scalar;
    mv
}

impl Div<&Real> for &CliffordNumber {
    type Output = CliffordNumber;
    fn div(self, rhs: &Real) -> CliffordNumber {
        let mut result = self.clone();
        for coefficient in result.coeffs_mut_slice() {
            *coefficient = &*coefficient / rhs;
        }
        result
    }
}
impl Div<&CliffordNumber> for &Real {
    type Output = CliffordNumber;
    fn div(self, rhs: &CliffordNumber) -> CliffordNumber {
        scale(&rhs.geometric_inverse(), self)
    }
}

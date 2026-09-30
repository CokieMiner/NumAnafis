//! [`Numeric`] trait implementation for [`CliffordNumber`].
//!
//! Forwards each algebraic operation to the spectral dispatcher,
//! mirroring the inherent kernels where closed forms exist.

use core::cmp::Ordering;

use crate::{
    traits::Numeric,
    types::{Complex, Number, Real},
};

use super::{BladeAlgebra, CliffordNumber, Cmplx, SpectralDispatch, SpectralFn};
// Structural ordering for deterministic sorting.

impl CliffordNumber {
    /// Compare with total ordering: generator sets first, then blade coefficients.
    ///
    /// Distinct algebras never compare `Equal`: generator-set lengths alone
    /// admit `A~B`, `B~C`, `A<C`, so length, IDs, and metrics lead.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        // Generator sets first (length, then IDs, then metrics), so distinct
        // algebras never compare Equal: lengths alone admit A~B, B~C, A<C.
        match self.gens.len().cmp(&other.gens.len()) {
            Ordering::Equal => {}
            ord @ (Ordering::Less | Ordering::Greater) => return ord,
        }
        for i in 0..self.gens.len() {
            match self.gens.id_at(i).cmp(&other.gens.id_at(i)) {
                Ordering::Equal => {}
                ord @ (Ordering::Less | Ordering::Greater) => return ord,
            }
            match self.gens.metric_at(i).cmp(&other.gens.metric_at(i)) {
                Ordering::Equal => {}
                ord @ (Ordering::Less | Ordering::Greater) => return ord,
            }
        }
        let limit = self.blade_count();
        for i in 0..limit {
            let cmp = self.coeff(i).total_cmp(other.coeff(i));
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        Ordering::Equal
    }

    /// Returns true if every blade coefficient is numerically zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        (0..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Returns true for the multiplicative identity (scalar one, rest zero).
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.coeff(0).is_one() && (1..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Returns true for negative one (scalar negative one, rest zero).
    #[must_use]
    pub fn is_neg_one(&self) -> bool {
        self.coeff(0).is_neg_one() && (1..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Whether this is a real integer with all non-scalar blades zero.
    #[must_use]
    pub fn is_int(&self) -> bool {
        self.coeff(0).is_int() && (1..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Whether every active blade coefficient is an integer.
    /// Finite integral floats pass; this does not establish computation exactness.
    #[must_use]
    pub fn is_int_ring(&self) -> bool {
        self.coeffs_slice().iter().all(Real::is_int)
    }

    /// Returns true for strictly negative scalars (zero non-scalar blades).
    #[must_use]
    pub fn is_negative(&self) -> bool {
        self.coeff(0).is_negative() && (1..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Returns true for strictly positive scalars (zero non-scalar blades).
    #[must_use]
    pub fn is_positive(&self) -> bool {
        self.coeff(0).is_positive() && (1..self.blade_count()).all(|i| self.coeff(i).is_zero())
    }

    /// Returns true if every blade coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        (0..self.blade_count()).all(|i| self.coeff(i).is_finite())
    }

    /// Projects every blade coefficient to the Float tier.
    #[must_use]
    pub fn to_float(&self) -> Self {
        let mut out = Self::zero_unchecked(self.gens.clone());
        for i in 0..self.blade_count() {
            out.coeffs_mut_slice()[i] = self.coeff(i).to_float();
        }
        out
    }

    /// Approximate equality using the maximum absolute coefficient distance.
    #[must_use]
    pub fn approx_eq_number(&self, other: &Self, tolerance: &Self) -> bool {
        if self.gens != other.gens
            || self.gens != tolerance.gens
            || !self.is_finite()
            || !other.is_finite()
            || !tolerance.is_finite()
        {
            return false;
        }
        let distance = self
            .coeffs_slice()
            .iter()
            .zip(other.coeffs_slice())
            .map(|(a, b)| (a - b).abs())
            .fold(
                Real::zero(),
                |largest, value| {
                    if value > largest { value } else { largest }
                },
            );
        let bound =
            tolerance
                .coeffs_slice()
                .iter()
                .map(Real::abs)
                .fold(
                    Real::zero(),
                    |largest, value| {
                        if value > largest { value } else { largest }
                    },
                );
        distance <= bound
    }

    /// Greater real scalar, retaining the selected operand's algebra context.
    /// Non-real or unordered values return NaN in this operand's algebra.
    #[must_use]
    pub fn max(&self, other: &Self) -> Self {
        match (self.central_scalar_parts(), other.central_scalar_parts()) {
            (Some((a, b)), Some((c, d))) if b.is_zero() && d.is_zero() => match a.partial_cmp(c) {
                Some(Ordering::Less) => other.clone(),
                Some(Ordering::Equal | Ordering::Greater) => self.clone(),
                None => Self::nan(self.gens.clone()),
            },
            _ => Self::nan(self.gens.clone()),
        }
    }

    /// Lesser real scalar, retaining the selected operand's algebra context.
    /// Non-real or unordered values return NaN in this operand's algebra.
    #[must_use]
    pub fn min(&self, other: &Self) -> Self {
        match (self.central_scalar_parts(), other.central_scalar_parts()) {
            (Some((a, b)), Some((c, d))) if b.is_zero() && d.is_zero() => match a.partial_cmp(c) {
                Some(Ordering::Greater) => other.clone(),
                Some(Ordering::Equal | Ordering::Less) => self.clone(),
                None => Self::nan(self.gens.clone()),
            },
            _ => Self::nan(self.gens.clone()),
        }
    }
}
// Numeric trait implementation for CliffordNumber

impl Numeric for CliffordNumber {
    fn sin(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sin)
    }
    fn cos(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Cos)
    }
    fn tan(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Tan)
    }
    fn cot(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Cot)
    }
    fn sec(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sec)
    }
    fn csc(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Csc)
    }

    fn asin(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Asin)
    }
    fn acos(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acos)
    }
    fn atan(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Atan)
    }
    fn acot(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acot)
    }
    fn asec(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Asec)
    }
    fn acsc(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acsc)
    }

    fn sinh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sinh)
    }
    fn cosh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Cosh)
    }
    fn tanh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Tanh)
    }
    fn coth(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Coth)
    }
    fn sech(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sech)
    }
    fn csch(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Csch)
    }

    fn asinh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Asinh)
    }
    fn acosh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acosh)
    }
    fn atanh(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Atanh)
    }
    fn acoth(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acoth)
    }
    fn asech(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Asech)
    }
    fn acsch(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Acsch)
    }

    fn exp(&self) -> Self {
        self.quadratic_exponential(false)
            .unwrap_or_else(|| SpectralDispatch::apply_via_spectral(self, &SpectralFn::Exp))
    }
    fn expm1(&self) -> Self {
        self.quadratic_exponential(true)
            .unwrap_or_else(|| SpectralDispatch::apply_via_spectral(self, &SpectralFn::Expm1))
    }
    fn exp_neg(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::ExpNeg)
    }
    fn ln(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Ln)
    }
    fn log1p(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Log1p)
    }

    fn sqrt(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sqrt)
    }
    fn cbrt(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Cbrt)
    }

    fn sinc(&self) -> Self {
        SpectralDispatch::apply_via_spectral(self, &SpectralFn::Sinc)
    }

    fn atan2(&self, x: &Self) -> Self {
        if let (Some((re, im)), Some((other_re, other_im))) =
            (self.central_scalar_parts(), x.central_scalar_parts())
        {
            let a = Number::from(Complex::from_parts(re.clone(), im.clone()));
            let b = Number::from(Complex::from_parts(other_re.clone(), other_im.clone()));
            return Self::scalar(self.gens.clone(), a.atan2(&b))
                .unwrap_or_else(|_| Self::nan(self.gens.clone()));
        }

        Self::nan(self.gens.clone())
    }

    #[cfg_attr(
        feature = "std",
        expect(
            clippy::suboptimal_flops,
            reason = "Clifford multivectors compute log_base via definition ln(self) / ln(base)"
        )
    )]
    fn log_base(&self, base: &Self) -> Self {
        if let (Some((re, im)), Some((other_re, other_im))) =
            (self.central_scalar_parts(), base.central_scalar_parts())
        {
            let a = Number::from(Complex::from_parts(re.clone(), im.clone()));
            let b = Number::from(Complex::from_parts(other_re.clone(), other_im.clone()));
            return Self::scalar(self.gens.clone(), a.log_base(&b))
                .unwrap_or_else(|_| Self::nan(self.gens.clone()));
        }

        self.ln() / base.ln()
    }

    fn pow(&self, exp: &Self) -> Self {
        if let (Some((re, im)), Some((other_re, other_im))) =
            (self.central_scalar_parts(), exp.central_scalar_parts())
        {
            let a = Number::from(Complex::from_parts(re.clone(), im.clone()));
            let b = Number::from(Complex::from_parts(other_re.clone(), other_im.clone()));
            return Self::scalar(self.gens.clone(), a.pow(&b))
                .unwrap_or_else(|_| Self::nan(self.gens.clone()));
        }

        if SpectralDispatch::is_pure_scalar(exp)
            && let Real::Int(exponent) = exp.coeff(0)
        {
            let mut factor = if *exponent < 0 {
                self.geometric_inverse()
            } else {
                self.clone()
            };
            let mut remaining = exponent.unsigned_abs();
            let mut result =
                Self::scalar(self.gens.clone(), Real::one()).expect("validated algebra");
            while remaining != 0 {
                if remaining & 1 != 0 {
                    result *= &factor;
                }
                remaining >>= 1;
                if remaining != 0 {
                    factor = &factor * &factor;
                }
            }
            return result;
        }
        if SpectralDispatch::is_pure_scalar(exp) {
            let exp_s = SpectralDispatch::extract_scalar(exp);
            let f = |c: &Cmplx| {
                let log_c = c.ln();
                log_c.scale(&exp_s).exp()
            };
            SpectralDispatch::apply_via_spectral_closure(self, &f)
        } else {
            (self.ln() * exp.clone()).exp()
        }
    }

    fn abs(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.abs())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn signum(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.signum())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn floor(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.floor())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn ceil(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.ceil())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn round(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.round())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn fract(&self) -> Self {
        if let Some((re, im)) = self.central_scalar_parts() {
            let scalar = Number::from(Complex::from_parts(re.clone(), im.clone()));
            Self::scalar(self.gens.clone(), scalar.fract())
                .unwrap_or_else(|_| Self::nan(self.gens.clone()))
        } else {
            Self::nan(self.gens.clone())
        }
    }
    fn negate(&self) -> Self {
        -self
    }
}

impl CliffordNumber {
    // For A=s+B with B²=q scalar, powers split into even and odd series.
    // exp(A)=exp(s)[C(q)+S(q)B]. This condition applies to any grade;
    // bivectors with non-scalar squares must retain the spectral path.
    fn quadratic_exponential(&self, minus_one: bool) -> Option<Self> {
        if !self.is_finite() {
            return None;
        }
        let scalar = self.coeff(0);
        let mut remainder = self.clone();
        *remainder.coeffs_mut_slice().first_mut()? = Real::zero();
        let mut active = remainder.nonzero_blades();
        let first = active.next();
        let (sign, magnitude) = match (first, active.next()) {
            (None, _) => (0, Real::zero()),
            (Some((blade, coefficient)), None) => {
                // A basis blade squares to its metric sign. Use |coefficient|
                // directly: squaring can overflow, underflow, or exceed rational
                // denominator capacity even when the exponential is representable.
                match BladeAlgebra::mul_blades(&self.gens, blade, blade) {
                    None => (0, Real::zero()),
                    Some((_, metric)) => (metric, coefficient.abs()),
                }
            }
            (Some(_), Some(_)) => {
                let square = remainder.geometric_mul(&remainder);
                if !SpectralDispatch::is_pure_scalar(&square) || !square.is_finite() {
                    return None;
                }
                let q = square.coeff(0);
                let metric = if q.is_zero() {
                    0
                } else if q.is_negative() {
                    -1
                } else {
                    1
                };
                (metric, q.abs().sqrt())
            }
        };
        drop(active);
        let (constant_factor, sine_ratio) = exponential_factors(sign, &magnitude, minus_one);
        let scale = scalar.exp();
        let constant = if minus_one {
            scalar.expm1() + &scale * &constant_factor
        } else {
            &scale * &constant_factor
        };
        let factor = scale * sine_ratio;
        if !constant.is_finite() || !factor.is_finite() {
            return None;
        }
        // Zero lanes stay zero, even when the approximate factor is large.
        for coefficient in remainder.coeffs_mut_slice() {
            if !coefficient.is_zero() {
                *coefficient = &*coefficient * &factor;
            }
        }
        *remainder.coeffs_mut_slice().first_mut()? = constant;
        Some(remainder)
    }
}

// Even/odd exponential series: elliptic (sign<0), hyperbolic (sign>0),
// and nilpotent (sign=0). Only expm1 needs the stable half-angle subtraction.
fn exponential_factors(sign: i8, magnitude: &Real, minus_one: bool) -> (Real, Real) {
    if sign == 0 {
        return (
            if minus_one { Real::zero() } else { Real::one() },
            Real::one(),
        );
    }
    let constant = if minus_one {
        // Approximate kernels do not need an exact half with a wider denominator.
        let half = magnitude.to_float() / Real::from_int(2);
        let half_sine = if sign < 0 { half.sin() } else { half.sinh() };
        let factor = if sign < 0 {
            Real::from_int(-2)
        } else {
            Real::from_int(2)
        };
        factor * &half_sine * half_sine
    } else if sign < 0 {
        magnitude.cos()
    } else {
        magnitude.cosh()
    };
    let linear = if sign < 0 {
        magnitude.sinc()
    } else {
        magnitude.sinh() / magnitude
    };
    (constant, linear)
}

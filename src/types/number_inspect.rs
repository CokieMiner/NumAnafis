//! Numeric inspection, explicit projections, and scalar component access.

use core::cmp::Ordering;

use crate::traits::Numeric;

use super::{Complex, FloatMath, IntType, Number, NumberRepr, Real};

macro_rules! scalar_predicate {
    ($name:ident, $test:expr, $doc:literal) => {
        #[doc = $doc]
        #[must_use]
        pub fn $name(&self) -> bool {
            self.scalar_parts().is_some_and($test)
        }
    };
}
impl Number {
    /// Whether the storage is a multivector.
    #[must_use]
    pub const fn is_clifford(&self) -> bool {
        match &self.0 {
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(_) => true,
            NumberRepr::Real(_) | NumberRepr::Complex(_) => false,
        }
    }
    /// Whether the rational representation is active.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        matches!(self.0, NumberRepr::Real(Real::Rational(_)))
    }
    /// Whether the real float representation is active.
    #[must_use]
    pub const fn is_float(&self) -> bool {
        matches!(self.0, NumberRepr::Real(Real::Float(_)))
    }
    scalar_predicate!(
        is_real,
        |(_, im)| im.is_zero(),
        "Whether this is an ordinary real scalar."
    );
    scalar_predicate!(
        is_complex,
        |(_, im)| !im.is_zero(),
        "Whether this is an ordinary scalar with nonzero imaginary part."
    );
    scalar_predicate!(
        is_int,
        |(re, im)| im.is_zero() && re.is_int(),
        "Whether this value is a real integer. Finite floats pass only when their stored value has no fractional part; this does not establish computation exactness."
    );
    scalar_predicate!(
        is_one,
        |(re, im)| im.is_zero() && re.is_one(),
        "Whether the value is the multiplicative identity."
    );
    scalar_predicate!(
        is_neg_one,
        |(re, im)| im.is_zero() && re.is_neg_one(),
        "Whether the value is negative one."
    );
    scalar_predicate!(
        is_negative,
        |(re, im)| im.is_zero() && re.is_negative(),
        "Whether the value is a strictly negative real scalar."
    );
    scalar_predicate!(
        is_positive,
        |(re, im)| im.is_zero() && re.is_positive(),
        "Whether the value is a strictly positive real scalar."
    );
    /// Whether every component is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        match &self.0 {
            NumberRepr::Real(v) => v.is_zero(),
            NumberRepr::Complex(v) => v.re.is_zero() && v.im.is_zero(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => v.is_zero(),
        }
    }
    /// Whether every real component is an integer, including integral floats.
    /// For complex scalars this tests membership in the Gaussian integers;
    /// for Clifford values it checks every active blade coefficient.
    /// This tests represented values, not the exactness of their computation.
    #[must_use]
    pub fn is_int_ring(&self) -> bool {
        !self.any_component(|component| !component.is_int())
    }
    /// Whether any component is NaN.
    #[must_use]
    pub fn is_nan(&self) -> bool {
        self.any_component(Real::is_nan)
    }
    /// Whether any component is infinite.
    #[must_use]
    pub fn is_infinite(&self) -> bool {
        self.any_component(Real::is_infinite)
    }
    /// Whether all components are finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        !self.is_nan() && !self.is_infinite()
    }
    /// Real component, or the grade-zero coefficient of a multivector.
    /// Other blades are excluded from this projection.
    #[must_use]
    pub fn re(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(value) => value.clone().into(),
            NumberRepr::Complex(value) => value.re.clone().into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(value) => value.coeff(0).clone().into(),
        }
    }
    /// Imaginary coefficient relative to the same unit used by ordinary complex
    /// scalars. For CGA this is I5. Other blades are excluded from this projection.
    /// Returns zero for real values or algebras without a declared complex unit.
    #[must_use]
    pub fn im(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(_) => Self::zero(),
            NumberRepr::Complex(value) => value.im.clone().into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(value) => value
                .generator_set()
                .complex_blade()
                .map_or_else(Self::zero, |blade| value.coeff(blade).clone().into()),
        }
    }
    /// Scalar conjugation or its compatible ambient anti-involution.
    /// Uses reversion when ordinary Clifford conjugation fixes the declared
    /// central imaginary unit; otherwise uses Clifford conjugation.
    #[must_use]
    pub fn conj(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(_) => self.clone(),
            NumberRepr::Complex(v) => v.conj().into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => {
                // In odd ranks 3 mod 4, grade conjugation fixes J. Reversion
                // instead sends J to -J and is a complex-compatible anti-involution.
                if v.generator_set().complex_blade().is_some() && v.n_generators() % 4 == 3 {
                    v.reverse().into()
                } else {
                    v.clifford_conjugate().into()
                }
            }
        }
    }
    /// Scalar squared modulus or the scalar product with the compatible
    /// ambient conjugate. General multivectors can have a signed form.
    #[must_use]
    pub fn norm_sq(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(v) => (v * v).into(),
            NumberRepr::Complex(v) => v.norm_sq().into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => v.central_scalar_parts().map_or_else(
                || {
                    if v.generator_set().complex_blade().is_some() && v.n_generators() % 4 == 3 {
                        v.scalar_product(&v.reverse()).into()
                    } else {
                        v.norm_sq().into()
                    }
                },
                |(re, im)| (re * re + im * im).into(),
            ),
        }
    }
    /// Principal scalar argument; returns NaN for a general multivector.
    #[must_use]
    pub fn arg(&self) -> Self {
        self.scalar_parts()
            .map_or_else(Self::nan, |(re, im)| im.atan2(re).into())
    }
    /// Explicit projection to the real scalar component as f64.
    ///
    /// # Panics
    /// Panics for a general multivector with no scalar complex interpretation.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        self.scalar_parts()
            .expect("scalar projection required")
            .0
            .to_f64()
    }
    /// Explicit real-component projection to f32.
    ///
    /// # Panics
    /// Panics for a general multivector.
    #[expect(
        clippy::as_conversions,
        clippy::cast_possible_truncation,
        reason = "Explicit narrowing float projection"
    )]
    #[must_use]
    pub fn to_f32(&self) -> f32 {
        self.to_f64() as f32
    }
    /// Truncating real-component projection to i64.
    ///
    /// # Panics
    /// Panics for a general multivector.
    #[must_use]
    pub fn to_int(&self) -> IntType {
        self.scalar_parts()
            .expect("scalar projection required")
            .0
            .to_int()
    }
    /// Explicit approximate conversion preserving every component.
    #[must_use]
    pub fn to_float(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(v) => v.to_float().into(),
            NumberRepr::Complex(v) => Complex::from_parts(v.re.to_float(), v.im.to_float()).into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => v.to_float().into(),
        }
    }
    /// Total order by scalar magnitude/components, then multivector context.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        match (self.scalar_parts(), other.scalar_parts()) {
            (Some((a, b)), Some((c, d))) => FloatMath::total_cmp(
                &FloatMath::hypot(&a.float_value(), &b.float_value()),
                &FloatMath::hypot(&c.float_value(), &d.float_value()),
            )
            .then_with(|| a.total_cmp(c))
            .then_with(|| b.total_cmp(d)),
            (Some(_), None) => Ordering::Less,
            (None, Some(_)) => Ordering::Greater,
            (None, None) => {
                #[cfg(feature = "clifford")]
                if let (NumberRepr::Clifford(a), NumberRepr::Clifford(b)) = (&self.0, &other.0) {
                    return a.total_cmp(b);
                }
                Ordering::Equal
            }
        }
    }
    /// Greater real scalar. Returns NaN for unordered or non-real operands.
    /// Use [`Self::total_cmp`] when an explicit sorting order is needed.
    #[must_use]
    pub fn max(&self, other: &Self) -> Self {
        if !self.is_real() || !other.is_real() {
            return Self::nan();
        }
        match self.partial_cmp(other) {
            Some(Ordering::Less) => other.clone(),
            Some(Ordering::Equal | Ordering::Greater) => self.clone(),
            None => Self::nan(),
        }
    }
    /// Lesser real scalar. Returns NaN for unordered or non-real operands.
    /// Use [`Self::total_cmp`] when an explicit sorting order is needed.
    #[must_use]
    pub fn min(&self, other: &Self) -> Self {
        if !self.is_real() || !other.is_real() {
            return Self::nan();
        }
        match self.partial_cmp(other) {
            Some(Ordering::Greater) => other.clone(),
            Some(Ordering::Equal | Ordering::Less) => self.clone(),
            None => Self::nan(),
        }
    }
    /// Applies the lower bound with max, then the upper bound with min.
    /// Non-real or unordered inputs return NaN. Bound ordering is not validated.
    #[must_use]
    pub fn clamp(&self, min: &Self, max: &Self) -> Self {
        self.max(min).min(max)
    }
    /// Finite scalar distance or matching-algebra coefficient distance.
    #[must_use]
    pub fn approx_eq_number(&self, other: &Self, tolerance: &Self) -> bool {
        #[cfg(feature = "clifford")]
        if let (NumberRepr::Clifford(a), NumberRepr::Clifford(b), NumberRepr::Clifford(t)) =
            (&self.0, &other.0, &tolerance.0)
        {
            return a.approx_eq_number(b, t);
        }
        if !self.is_finite() || !other.is_finite() || !tolerance.is_finite() {
            return false;
        }
        match ((self - other).scalar_parts(), tolerance.scalar_parts()) {
            (Some((re, im)), Some((tr, ti))) => {
                FloatMath::hypot(&re.float_value(), &im.float_value())
                    <= FloatMath::hypot(&tr.float_value(), &ti.float_value())
            }
            _ => false,
        }
    }
    fn any_component(&self, predicate: fn(&Real) -> bool) -> bool {
        match &self.0 {
            NumberRepr::Real(v) => predicate(v),
            NumberRepr::Complex(v) => predicate(&v.re) || predicate(&v.im),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => v.coeffs_slice().iter().any(predicate),
        }
    }
}

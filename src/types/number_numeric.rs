//! Domain promotion and elementary evaluator dispatch.

#[cfg(feature = "clifford")]
use crate::clifford::CliffordNumber;
use crate::traits::Numeric;

use super::{Complex, Number, NumberRepr, Real};

macro_rules! unary {
    ($($name:ident),*) => { $(fn $name(&self)->Self {
        match &self.0 {
            NumberRepr::Real(value)=>value.$name().into(),
            NumberRepr::Complex(value)=>Numeric::$name(value.as_ref()).into(),
            #[cfg(feature="clifford")] NumberRepr::Clifford(value)=>value.$name().into(),
        }
    })* };
}
macro_rules! domain_unary {
    ($name:ident,$outside:expr) => {
        fn $name(&self) -> Self {
            match &self.0 {
                NumberRepr::Real(value) if ($outside)(value) => {
                    Numeric::$name(&Complex::from_parts(value.clone(), Real::zero())).into()
                }
                NumberRepr::Real(value) => value.$name().into(),
                NumberRepr::Complex(value) => Numeric::$name(value.as_ref()).into(),
                #[cfg(feature = "clifford")]
                NumberRepr::Clifford(value) => value.$name().into(),
            }
        }
    };
}
macro_rules! binary {
    ($name:ident) => {
        fn $name(&self, other: &Self) -> Self {
            #[cfg(feature = "clifford")]
            if self.is_clifford() || other.is_clifford() {
                let context = self
                    .as_clifford()
                    .or_else(|| other.as_clifford())
                    .expect("Clifford branch has a context");
                let a = match &self.0 {
                    NumberRepr::Clifford(v) => (**v).clone(),
                    NumberRepr::Real(_) | NumberRepr::Complex(_) => {
                        match CliffordNumber::scalar(context.generator_set().clone(), self.clone())
                        {
                            Ok(v) => v,
                            Err(_) => {
                                return CliffordNumber::nan(context.generator_set().clone()).into();
                            }
                        }
                    }
                };
                let b = match &other.0 {
                    NumberRepr::Clifford(v) => (**v).clone(),
                    NumberRepr::Real(_) | NumberRepr::Complex(_) => {
                        match CliffordNumber::scalar(context.generator_set().clone(), other.clone())
                        {
                            Ok(v) => v,
                            Err(_) => {
                                return CliffordNumber::nan(context.generator_set().clone()).into();
                            }
                        }
                    }
                };
                return a.$name(&b).into();
            }
            self.scalar_binary(other, |a, b| a.$name(b), |a, b| Numeric::$name(a, b))
        }
    };
}
impl Numeric for Number {
    unary!(
        sin, cos, tan, sinh, cosh, tanh, exp, expm1, exp_neg, cbrt, signum, floor, ceil, round,
        fract, negate, sinc, atan, asinh
    );
    domain_unary!(sqrt, Real::is_negative);
    domain_unary!(ln, Real::is_negative);
    domain_unary!(log1p, |v: &Real| v < &Real::neg_one());
    domain_unary!(asin, |v: &Real| v.abs() > Real::one());
    domain_unary!(acos, |v: &Real| v.abs() > Real::one());
    domain_unary!(acosh, |v: &Real| v < &Real::one());
    domain_unary!(atanh, |v: &Real| v.abs() > Real::one());
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
    fn abs(&self) -> Self {
        match &self.0 {
            NumberRepr::Real(v) => v.abs().into(),
            NumberRepr::Complex(v) => v.abs().into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => v.abs().into(),
        }
    }
    binary!(atan2);
    binary!(log_base);
    fn pow(&self, exponent: &Self) -> Self {
        if let NumberRepr::Real(Real::Rational(value)) = &exponent.0
            && *value.denom() == 2
        {
            return self.sqrt().pow(&Self::from_int(*value.numer()));
        }

        #[cfg(feature = "clifford")]
        if self.is_clifford() || exponent.is_clifford() {
            return self.clifford_power(exponent);
        }
        if self.is_negative() && !exponent.is_int() {
            return self
                .clone()
                .into_complex()
                .expect("scalar power")
                .pow(&exponent.clone().into_complex().expect("scalar exponent"))
                .into();
        }
        self.scalar_binary(exponent, Real::pow, Complex::pow)
    }
}
impl Number {
    /// Truncates the scalar components; unsupported Clifford inputs return NaN.
    #[must_use]
    pub fn trunc(&self) -> Self {
        self.clone() - self.fract()
    }
    fn scalar_binary(
        &self,
        other: &Self,
        real: fn(&Real, &Real) -> Real,
        complex: fn(&Complex, &Complex) -> Complex,
    ) -> Self {
        if let (NumberRepr::Real(a), NumberRepr::Real(b)) = (&self.0, &other.0) {
            return real(a, b).into();
        }
        self.with_complex(|a| other.with_complex(|b| complex(a, b).into()))
    }
    #[cfg(feature = "clifford")]
    fn clifford_power(&self, exponent: &Self) -> Self {
        let context = self
            .as_clifford()
            .or_else(|| exponent.as_clifford())
            .expect("Clifford power context");
        let a = match &self.0 {
            NumberRepr::Clifford(v) => (**v).clone(),
            NumberRepr::Real(_) | NumberRepr::Complex(_) => {
                match CliffordNumber::scalar(context.generator_set().clone(), self.clone()) {
                    Ok(v) => v,
                    Err(_) => return CliffordNumber::nan(context.generator_set().clone()).into(),
                }
            }
        };
        let b = match &exponent.0 {
            NumberRepr::Clifford(v) => (**v).clone(),
            NumberRepr::Real(_) | NumberRepr::Complex(_) => {
                match CliffordNumber::scalar(context.generator_set().clone(), exponent.clone()) {
                    Ok(v) => v,
                    Err(_) => return CliffordNumber::nan(context.generator_set().clone()).into(),
                }
            }
        };
        a.pow(&b).into()
    }
}

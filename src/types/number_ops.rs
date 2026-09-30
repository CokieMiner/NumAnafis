//! Direct tier arithmetic and exact value semantics.

use core::{
    cmp::Ordering,
    fmt::{Debug, Display, Formatter, Result as FmtResult},
    hash::{Hash, Hasher},
    ops::{Add, Div, Mul, Neg, Sub},
};

#[cfg(feature = "clifford")]
use crate::clifford::CliffordNumber;

use super::{Complex, Number, NumberRepr};

macro_rules! binary_arithmetic {
    ($trait:ident,$method:ident,$operator:tt) => {
        impl $trait for &Number {
            type Output=Number;
            fn $method(self,rhs:Self)->Number {
                if let (NumberRepr::Real(a),NumberRepr::Real(b))=(&self.0,&rhs.0) { return (a $operator b).into(); }
                #[cfg(feature="clifford")]
                if self.is_clifford() || rhs.is_clifford() {
                    return match (&self.0,&rhs.0) {
                        (NumberRepr::Clifford(a),NumberRepr::Clifford(b))=>{
                            if a.generator_set()!=b.generator_set() {
                                if let Ok(scalar)=self.extract_complex() { return Number::from(scalar).$method(rhs); }
                                if let Ok(scalar)=rhs.extract_complex() { return self.$method(&Number::from(scalar)); }
                            }
                            (a.as_ref() $operator b.as_ref()).into()
                        },
                        (NumberRepr::Clifford(a),NumberRepr::Real(_) | NumberRepr::Complex(_))=>match CliffordNumber::scalar(a.generator_set().clone(),rhs.clone()) {
                            Ok(b)=>(a.as_ref() $operator &b).into(),Err(_)=>CliffordNumber::nan(a.generator_set().clone()).into(),
                        },
                        (NumberRepr::Real(_) | NumberRepr::Complex(_),NumberRepr::Clifford(b))=>match CliffordNumber::scalar(b.generator_set().clone(),self.clone()) {
                            Ok(a)=>(&a $operator b.as_ref()).into(),Err(_)=>CliffordNumber::nan(b.generator_set().clone()).into(),
                        },
                        (NumberRepr::Real(_) | NumberRepr::Complex(_),NumberRepr::Real(_) | NumberRepr::Complex(_))=>Number::nan(),
                    };
                }
                self.with_complex(|a|rhs.with_complex(|b|Complex::$method(a,b).into()))
            }
        }
        impl $trait for Number { type Output=Self; fn $method(self,rhs:Self)->Self { (&self).$method(&rhs) } }
        impl $trait<&Self> for Number { type Output=Self; fn $method(self,rhs:&Self)->Self { (&self).$method(rhs) } }
        impl $trait<Number> for &Number { type Output=Number; fn $method(self,rhs:Number)->Number { self.$method(&rhs) } }
    };
}
binary_arithmetic!(Add,add,+);
binary_arithmetic!(Sub,sub,-);
binary_arithmetic!(Mul,mul,*);
binary_arithmetic!(Div,div,/);

// Primitive literals use the same dispatch as Number operands. Keep operand
// order explicit so subtraction/division and Clifford products retain meaning.
macro_rules! literal_arithmetic {
    ($trait:ident, $method:ident, $operator:tt; $($ty:ty),*) => {
        $(
            impl $trait<$ty> for Number {
                type Output = Self;
                fn $method(self, rhs: $ty) -> Self { &self $operator &Self::from(rhs) }
            }
            impl $trait<$ty> for &Number {
                type Output = Number;
                fn $method(self, rhs: $ty) -> Number { self $operator &Number::from(rhs) }
            }
            impl $trait<Number> for $ty {
                type Output = Number;
                fn $method(self, rhs: Number) -> Number { &Number::from(self) $operator &rhs }
            }
            impl $trait<&Number> for $ty {
                type Output = Number;
                fn $method(self, rhs: &Number) -> Number { &Number::from(self) $operator rhs }
            }
        )*
    };
}
literal_arithmetic!(Add, add, +; i64, f64);
literal_arithmetic!(Sub, sub, -; i64, f64);
literal_arithmetic!(Mul, mul, *; i64, f64);
literal_arithmetic!(Div, div, /; i64, f64);

impl Neg for Number {
    type Output = Self;
    fn neg(self) -> Self {
        match self.0 {
            NumberRepr::Real(v) => (-v).into(),
            NumberRepr::Complex(v) => (-*v).into(),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => (-*v).into(),
        }
    }
}
impl Neg for &Number {
    type Output = Number;
    fn neg(self) -> Number {
        -self.clone()
    }
}
impl PartialEq for Number {
    fn eq(&self, other: &Self) -> bool {
        if let (Some((a, b)), Some((c, d))) = (self.scalar_parts(), other.scalar_parts()) {
            return a == c && b == d;
        }
        #[cfg(feature = "clifford")]
        if let (NumberRepr::Clifford(a), NumberRepr::Clifford(b)) = (&self.0, &other.0) {
            return a == b;
        }
        false
    }
}
impl PartialOrd for Number {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if let (Some((a, b)), Some((c, d))) = (self.scalar_parts(), other.scalar_parts()) {
            return if b.is_zero() && d.is_zero() {
                a.partial_cmp(c)
            } else {
                (self == other).then_some(Ordering::Equal)
            };
        }
        #[cfg(feature = "clifford")]
        if let (NumberRepr::Clifford(a), NumberRepr::Clifford(b)) = (&self.0, &other.0) {
            return a.partial_cmp(b);
        }
        None
    }
}
impl Hash for Number {
    fn hash<H: Hasher>(&self, state: &mut H) {
        if let Some((re, im)) = self.scalar_parts() {
            0_u8.hash(state);
            re.hash(state);
            im.hash(state);
        } else {
            #[cfg(feature = "clifford")]
            if let NumberRepr::Clifford(value) = &self.0 {
                1_u8.hash(state);
                value.hash(state);
            }
        }
    }
}
impl Display for Number {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match &self.0 {
            NumberRepr::Real(v) => Display::fmt(v, f),
            NumberRepr::Complex(v) => Display::fmt(v, f),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => Display::fmt(v, f),
        }
    }
}
impl Debug for Number {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match &self.0 {
            NumberRepr::Real(v) => Debug::fmt(v, f),
            NumberRepr::Complex(v) => Debug::fmt(v, f),
            #[cfg(feature = "clifford")]
            NumberRepr::Clifford(v) => Debug::fmt(v, f),
        }
    }
}

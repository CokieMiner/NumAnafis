//! Spectral function AST for Clifford functional calculus.
//!
//! The [`SpectralFn`] vocabulary (trigonometric, hyperbolic, exponential,
//! and composite nodes) with eigenvalue evaluation and
//! symbolic differentiation. Series evaluation lives in
//! [`super::spectral_jet`].

use alloc::{boxed::Box, vec, vec::Vec};

use crate::types::{RationalType, Real};

use super::Cmplx;

#[derive(Debug, Clone, PartialEq)]
pub enum SpectralFn {
    // Trigonometric
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Acot,
    Asec,
    Acsc,

    // Hyperbolic
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Sech,
    Csch,
    Asinh,
    Acosh,
    Atanh,
    Acoth,
    Asech,
    Acsch,

    // Exponential & Power
    Exp,
    Expm1,
    ExpNeg,
    Ln,
    Log1p,
    Sqrt,
    Cbrt,

    // Cardinal sine
    Sinc,

    // AST Nodes
    Var,
    Const(Cmplx),
    Add(Vec<Self>),
    Mul(Vec<Self>),
    Inv(Box<Self>),
    Neg(Box<Self>),
    Compose(Box<Self>, Box<Self>),
}

impl SpectralFn {
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "match ergonomics produce cleaner code here"
    )]
    pub fn eval(&self, c: &Cmplx) -> Cmplx {
        match self {
            Self::Sin => c.sin(),
            Self::Cos => c.cos(),
            Self::Tan => c.sin().div(&c.cos()),
            Self::Cot => c.cos().div(&c.sin()),
            Self::Sec => Cmplx::one().div(&c.cos()),
            Self::Csc => Cmplx::one().div(&c.sin()),
            Self::Asin => c.asin(),
            Self::Acos => c.acos(),
            Self::Atan => c.atan(),
            Self::Acot => c.acot(),
            Self::Asec => c.asec(),
            Self::Acsc => c.acsc(),

            Self::Sinh => c.sinh(),
            Self::Cosh => c.cosh(),
            Self::Tanh => c.sinh().div(&c.cosh()),
            Self::Coth => c.cosh().div(&c.sinh()),
            Self::Sech => Cmplx::one().div(&c.cosh()),
            Self::Csch => Cmplx::one().div(&c.sinh()),
            Self::Asinh => c.asinh(),
            Self::Acosh => c.acosh(),
            Self::Atanh => c.atanh(),
            Self::Acoth => c.acoth(),
            Self::Asech => c.asech(),
            Self::Acsch => c.acsch(),

            Self::Exp => c.exp(),
            Self::Expm1 => c.expm1(),
            Self::ExpNeg => c.neg().exp(),
            Self::Ln => c.ln(),
            Self::Log1p => c.log1p(),
            Self::Sqrt => c.sqrt(),
            Self::Cbrt => c.cbrt(),

            Self::Sinc => {
                if c.0.is_zero() && c.1.is_zero() {
                    Cmplx::one()
                } else {
                    c.sin().div(c)
                }
            }

            Self::Var => c.clone(),
            Self::Const(v) => v.clone(),
            Self::Add(args) => {
                let mut sum = Cmplx::zero();
                for a in args {
                    sum = sum.add(&a.eval(c));
                }
                sum
            }
            Self::Mul(args) => {
                let mut prod = Cmplx::one();
                for a in args {
                    prod = prod.mul(&a.eval(c));
                }
                prod
            }
            Self::Inv(a) => Cmplx::one().div(&a.eval(c)),
            Self::Neg(a) => a.eval(c).neg(),
            Self::Compose(outer, inner) => outer.eval(&inner.eval(c)),
        }
    }

    #[expect(
        clippy::too_many_lines,
        clippy::pattern_type_mismatch,
        reason = "Large match block and match ergonomics are much cleaner here"
    )]
    pub fn derivative(&self) -> Self {
        match self {
            Self::Sin => Self::Cos,
            Self::Cos => Self::Neg(Box::new(Self::Sin)),
            Self::Tan => Self::Mul(vec![Self::Sec, Self::Sec]),
            Self::Cot => Self::Neg(Box::new(Self::Mul(vec![Self::Csc, Self::Csc]))),
            Self::Sec => Self::Mul(vec![Self::Sec, Self::Tan]),
            Self::Csc => Self::Neg(Box::new(Self::Mul(vec![Self::Csc, Self::Cot]))),

            Self::Asin => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Const(Cmplx::one()),
                    Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                ])),
            ))),
            Self::Acos => Self::Neg(Box::new(Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Const(Cmplx::one()),
                    Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                ])),
            ))))),
            Self::Atan => Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Mul(vec![Self::Var, Self::Var]),
            ]))),
            Self::Acot => Self::Neg(Box::new(Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Mul(vec![Self::Var, Self::Var]),
            ]))))),
            Self::Asec => Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Mul(vec![Self::Var, Self::Var]),
                        Self::Const(Cmplx::one().neg()),
                    ])),
                ),
            ]))),
            Self::Acsc => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Mul(vec![Self::Var, Self::Var]),
                        Self::Const(Cmplx::one().neg()),
                    ])),
                ),
            ]))))),

            Self::Sinh => Self::Cosh,
            Self::Cosh => Self::Sinh,
            Self::Tanh => Self::Mul(vec![Self::Sech, Self::Sech]),
            Self::Coth => Self::Neg(Box::new(Self::Mul(vec![Self::Csch, Self::Csch]))),
            Self::Sech => Self::Neg(Box::new(Self::Mul(vec![Self::Sech, Self::Tanh]))),
            Self::Csch => Self::Neg(Box::new(Self::Mul(vec![Self::Csch, Self::Coth]))),

            Self::Asinh => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Mul(vec![Self::Var, Self::Var]),
                    Self::Const(Cmplx::one()),
                ])),
            ))),
            Self::Acosh => Self::Inv(Box::new(Self::Compose(
                Box::new(Self::Sqrt),
                Box::new(Self::Add(vec![
                    Self::Mul(vec![Self::Var, Self::Var]),
                    Self::Const(Cmplx::one().neg()),
                ])),
            ))),
            Self::Atanh | Self::Acoth => Self::Inv(Box::new(Self::Add(vec![
                Self::Const(Cmplx::one()),
                Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
            ]))),
            Self::Asech => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Const(Cmplx::one()),
                        Self::Neg(Box::new(Self::Mul(vec![Self::Var, Self::Var]))),
                    ])),
                ),
            ]))))),
            Self::Acsch => Self::Neg(Box::new(Self::Inv(Box::new(Self::Mul(vec![
                Self::Var,
                Self::Compose(
                    Box::new(Self::Sqrt),
                    Box::new(Self::Add(vec![
                        Self::Const(Cmplx::one()),
                        Self::Mul(vec![Self::Var, Self::Var]),
                    ])),
                ),
            ]))))),

            Self::Exp | Self::Expm1 => Self::Exp,
            Self::ExpNeg => Self::Neg(Box::new(Self::ExpNeg)),
            Self::Ln => Self::Inv(Box::new(Self::Var)),
            Self::Log1p => Self::Inv(Box::new(Self::Add(vec![
                Self::Var,
                Self::Const(Cmplx::one()),
            ]))),

            Self::Sqrt => Self::Mul(vec![
                Self::Const(Cmplx::new(
                    Real::from_rational(RationalType::from_parts(1, 2)),
                    Real::zero(),
                )),
                Self::Inv(Box::new(Self::Sqrt)),
            ]),
            Self::Cbrt => Self::Mul(vec![
                Self::Const(Cmplx::new(
                    Real::from_rational(RationalType::from_parts(1, 3)),
                    Real::zero(),
                )),
                Self::Inv(Box::new(Self::Mul(vec![Self::Cbrt, Self::Cbrt]))),
            ]),

            Self::Sinc => Self::Mul(vec![
                Self::Inv(Box::new(Self::Var)),
                Self::Add(vec![Self::Cos, Self::Neg(Box::new(Self::Sinc))]),
            ]),

            // Base cases
            Self::Var => Self::Const(Cmplx::one()),
            Self::Const(_) => Self::Const(Cmplx::zero()),

            // AST Rules
            Self::Compose(outer, inner) => Self::Mul(vec![
                Self::Compose(Box::new(outer.derivative()), inner.clone()),
                inner.derivative(),
            ]),
            Self::Neg(a) => Self::Neg(Box::new(a.derivative())),
            Self::Inv(a) => Self::Neg(Box::new(Self::Mul(vec![
                a.derivative(),
                Self::Inv(Box::new(Self::Mul(vec![*a.clone(), *a.clone()]))),
            ]))),
            Self::Add(args) => Self::Add(args.iter().map(Self::derivative).collect()),
            Self::Mul(args) => {
                let mut sum_terms = Vec::with_capacity(args.len());
                for i in 0..args.len() {
                    let mut prod_terms = args.clone();
                    prod_terms[i] = prod_terms[i].derivative();
                    sum_terms.push(Self::Mul(prod_terms));
                }
                Self::Add(sum_terms)
            }
        }
    }
}

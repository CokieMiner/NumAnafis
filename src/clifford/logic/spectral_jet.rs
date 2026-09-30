//! Jet (truncated Taylor series) evaluation of [`SpectralFn`].
//!
//! Dual-number Taylor arithmetic plus the [`SpectralFn::eval_jet`]
//! interpreter, which evaluates an AST once over a series instead of
//! differentiating in a loop.

use alloc::{vec, vec::Vec};

use crate::{
    traits::Numeric,
    types::{IntType, Real},
};

use super::{Cmplx, SpectralFn};
// Jet — truncated Taylor series (dual numbers of degree K)
//
// Each Jet represents a truncated Taylor series
//   f(x₀ + ε) = f(x₀) + f'(x₀) ε + f''(x₀)/2 ε² + … + f^(K)(x₀)/K! εᴷ
// where coeffs[i] = f^(i)(x₀) / i! .
//
// Evaluated once instead of calling `.derivative()` in a loop (which can
// grow the AST exponentially via the product/chain rules).

#[derive(Debug, Clone, PartialEq)]
pub struct Jet {
    pub coeffs: Vec<Cmplx>,
}

impl Jet {
    pub const fn new(coeffs: Vec<Cmplx>) -> Self {
        Self { coeffs }
    }

    pub const fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    /// Constant jet of the given degree.
    ///
    /// # Panics
    /// Panics if `degree` is `usize::MAX` (the coefficient count overflows).
    pub fn constant(c: Cmplx, degree: usize) -> Self {
        let len = degree.checked_add(1).expect("jet degree overflows usize");
        let mut coeffs = vec![Cmplx::zero(); len];
        coeffs[0] = c;
        Self { coeffs }
    }

    pub fn zero(degree: usize) -> Self {
        Self::constant(Cmplx::zero(), degree)
    }

    pub fn one(degree: usize) -> Self {
        Self::constant(Cmplx::one(), degree)
    }

    /// Variable jet: `λ + ε` (used as the differentiation point in Block Parlett).
    ///
    /// # Panics
    /// Panics if `degree` is `usize::MAX` (the coefficient count overflows).
    pub fn variable(lambda: &Cmplx, degree: usize) -> Self {
        let len = degree.checked_add(1).expect("jet degree overflows usize");
        let mut coeffs = vec![Cmplx::zero(); len];
        coeffs[0] = lambda.clone();
        if degree >= 1 {
            coeffs[1] = Cmplx::one();
        }
        Self { coeffs }
    }
    // Arithmetic

    pub fn add(&self, other: &Self) -> Self {
        let max_deg = self.degree().max(other.degree());
        let mut coeffs = vec![Cmplx::zero(); max_deg + 1];
        for (i, c) in coeffs.iter_mut().enumerate().take(self.degree() + 1) {
            *c = c.add(&self.coeffs[i]);
        }
        for (i, c) in coeffs.iter_mut().enumerate().take(other.degree() + 1) {
            *c = c.add(&other.coeffs[i]);
        }
        Self { coeffs }
    }

    pub fn neg(&self) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(Cmplx::neg).collect(),
        }
    }

    /// Scale every coefficient by a complex scalar.
    pub fn scale(&self, scalar: &Cmplx) -> Self {
        Self {
            coeffs: self.coeffs.iter().map(|c| c.mul(scalar)).collect(),
        }
    }

    /// Cauchy product of two truncated Taylor series:
    ///   (f·g)[n] = Σ_{i+j=n} f[i]·g[j]
    ///
    /// Computation is modulo epsilon^(K+1), K the larger requested degree.
    /// Higher terms cannot contribute to a coefficient of order at most K.
    pub fn mul(&self, other: &Self) -> Self {
        let k = self.degree().max(other.degree());
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        for i in 0..=self.degree() {
            for j in 0..=other.degree().min(k - i) {
                coeffs[i + j] += &self.coeffs[i].mul(&other.coeffs[j]);
            }
        }
        Self { coeffs }
    }

    /// Multiplicative inverse via recurrence:
    ///   inv[0] = 1 / f[0]
    ///   inv[n] = -(1/f[0]) · Σ_{i=1}^{n} f[i] · inv[n-i]
    pub fn inv(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        let inv0 = Cmplx::one().div(&self.coeffs[0]);
        coeffs[0] = inv0.clone();
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..=n {
                sum += &self.coeffs[i].mul(&coeffs[n - i]);
            }
            coeffs[n] = sum.neg().mul(&inv0);
        }
        Self { coeffs }
    }

    /// Division: `a / b = a · b⁻¹`
    pub fn div(&self, other: &Self) -> Self {
        self.mul(&other.inv())
    }
    // Transcendental recurrences

    /// Exponential via recurrence:
    ///   exp[0] = exp(f[0])
    ///   exp[n] = (1/n) · Σ_{i=1}^{n} i · f[i] · exp[n-i]
    pub fn exp(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].exp();
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(
                    Real::from_int(
                        IntType::try_from(i).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                );
                sum += &i_c.mul(&self.coeffs[i]).mul(&coeffs[n - i]);
            }
            let n_c = Cmplx::new(
                Real::from_int(IntType::try_from(n).expect("Jet degree fits the integer backend")),
                Real::zero(),
            );
            coeffs[n] = sum.div(&n_c);
        }
        Self { coeffs }
    }

    /// Natural log via recurrence:
    ///   ln[0] = ln(f[0])
    ///   ln[n] = (f[n] - (1/n) · Σ_{i=1}^{n-1} i · ln[i] · f[n-i]) / f[0]
    pub fn ln(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].ln();
        let inv_f0 = Cmplx::one().div(&self.coeffs[0]);
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for (i, ci) in coeffs.iter().enumerate().take(n).skip(1) {
                let i_c = Cmplx::new(
                    Real::from_int(
                        IntType::try_from(i).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                );
                sum += &i_c.mul(ci).mul(&self.coeffs[n - i]);
            }
            coeffs[n] = self.coeffs[n]
                .sub(&sum.div(&Cmplx::new(
                    Real::from_int(
                        IntType::try_from(n).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                )))
                .mul(&inv_f0);
        }
        Self { coeffs }
    }

    /// Sine (and simultaneously cosine) via mutual recurrence:
    ///   sin[0] = sin(f[0])   cos[0] = cos(f[0])
    ///   sin[n] = (1/n) · Σ i·f[i]·cos[n-i]
    ///   cos[n] = -(1/n) · Σ i·f[i]·sin[n-i]
    fn sin_cos(&self) -> (Self, Self) {
        let k = self.degree();
        let mut sin_c = vec![Cmplx::zero(); k + 1];
        let mut cos_c = vec![Cmplx::zero(); k + 1];
        sin_c[0] = self.coeffs[0].sin();
        cos_c[0] = self.coeffs[0].cos();
        for n in 1..=k {
            let mut s_sum = Cmplx::zero();
            let mut c_sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(
                    Real::from_int(
                        IntType::try_from(i).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                );
                let fi = self.coeffs[i].mul(&i_c);
                s_sum = s_sum.add(&fi.mul(&cos_c[n - i]));
                c_sum = c_sum.add(&fi.mul(&sin_c[n - i]));
            }
            let n_c = Cmplx::new(
                Real::from_int(IntType::try_from(n).expect("Jet degree fits the integer backend")),
                Real::zero(),
            );
            sin_c[n] = s_sum.div(&n_c);
            cos_c[n] = c_sum.neg().div(&n_c);
        }
        (Self { coeffs: sin_c }, Self { coeffs: cos_c })
    }

    pub fn sin(&self) -> Self {
        self.sin_cos().0
    }

    pub fn cos(&self) -> Self {
        self.sin_cos().1
    }

    /// Hyperbolic sine (and simultaneously hyperbolic cosine) via mutual recurrence:
    ///   sinh[0] = sinh(f[0])  cosh[0] = cosh(f[0])
    ///   sinh[n] = (1/n) · Σ i·f[i]·cosh[n-i]
    ///   cosh[n] = (1/n) · Σ i·f[i]·sinh[n-i]
    fn sinh_cosh(&self) -> (Self, Self) {
        let k = self.degree();
        let mut sinh_c = vec![Cmplx::zero(); k + 1];
        let mut cosh_c = vec![Cmplx::zero(); k + 1];
        sinh_c[0] = self.coeffs[0].sinh();
        cosh_c[0] = self.coeffs[0].cosh();
        for n in 1..=k {
            let mut s_sum = Cmplx::zero();
            let mut c_sum = Cmplx::zero();
            for i in 1..=n {
                let i_c = Cmplx::new(
                    Real::from_int(
                        IntType::try_from(i).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                );
                let fi = self.coeffs[i].mul(&i_c);
                s_sum = s_sum.add(&fi.mul(&cosh_c[n - i]));
                c_sum = c_sum.add(&fi.mul(&sinh_c[n - i]));
            }
            let n_c = Cmplx::new(
                Real::from_int(IntType::try_from(n).expect("Jet degree fits the integer backend")),
                Real::zero(),
            );
            sinh_c[n] = s_sum.div(&n_c);
            cosh_c[n] = c_sum.div(&n_c);
        }
        (Self { coeffs: sinh_c }, Self { coeffs: cosh_c })
    }

    pub fn sinh(&self) -> Self {
        self.sinh_cosh().0
    }

    pub fn cosh(&self) -> Self {
        self.sinh_cosh().1
    }

    /// Square root via recurrence:
    ///   sqrt[0] = sqrt(f[0])
    ///   sqrt[n] = (f[n] - Σ_{i=1}^{n-1} sqrt[i]·sqrt[n-i]) / (2·sqrt[0])
    pub fn sqrt(&self) -> Self {
        let k = self.degree();
        let mut coeffs = vec![Cmplx::zero(); k + 1];
        coeffs[0] = self.coeffs[0].sqrt();
        let two_s0 = coeffs[0].add(&coeffs[0]); // 2·sqrt(f0)
        for n in 1..=k {
            let mut sum = Cmplx::zero();
            for i in 1..n {
                sum += &coeffs[i].mul(&coeffs[n - i]);
            }
            coeffs[n] = self.coeffs[n].sub(&sum).div(&two_s0);
        }
        Self { coeffs }
    }
}

impl SpectralFn {
    /// Evaluate this AST over a **Jet** (truncated Taylor series).
    ///
    /// The result is a Jet whose coefficients are the Taylor coefficients of the
    /// function composed with the input series, obviating repeated calls to
    /// `.derivative()`.
    #[expect(
        clippy::pattern_type_mismatch,
        clippy::too_many_lines,
        reason = "match ergonomics produce cleaner code here, and a centralized match block is more readable than fragmented helper functions"
    )]
    pub fn eval_jet(&self, x: &Jet) -> Jet {
        let deg = x.degree();
        match self {
            // --- Primitives with direct Jet recurrences ---
            Self::Exp => x.exp(),
            Self::Expm1 => {
                let mut result = x.exp();
                result.coeffs[0] = x.coeffs[0].expm1();
                result
            }
            Self::ExpNeg => x.neg().exp(),
            Self::Ln => x.ln(),
            Self::Log1p => {
                let one = Jet::one(deg);
                let mut result = x.add(&one).ln();
                result.coeffs[0] = x.coeffs[0].log1p();
                result
            }
            Self::Sin => x.sin(),
            Self::Cos => x.cos(),
            Self::Tan => {
                let scalar = x.sin();
                let c = x.cos();
                scalar.div(&c)
            }
            Self::Cot => {
                let scalar = x.sin();
                let c = x.cos();
                c.div(&scalar)
            }
            Self::Sec => Jet::one(deg).div(&x.cos()),
            Self::Csc => Jet::one(deg).div(&x.sin()),
            Self::Sinh => x.sinh(),
            Self::Cosh => x.cosh(),
            Self::Tanh => {
                let scalar = x.sinh();
                let c = x.cosh();
                scalar.div(&c)
            }
            Self::Coth => {
                let scalar = x.sinh();
                let c = x.cosh();
                c.div(&scalar)
            }
            Self::Sech => Jet::one(deg).div(&x.cosh()),
            Self::Csch => Jet::one(deg).div(&x.sinh()),

            Self::Sqrt => x.sqrt(),
            Self::Cbrt => {
                if x.coeffs[0].1.is_zero() {
                    // Real branch (matches scalar Real::cbrt): y₀ is the
                    // real cube root, higher coefficients from y³ = x.
                    Self::cbrt_real_jet(x)
                } else {
                    // Non-real inputs use exp(ln(x) / 3) with an exact rational exponent.
                    let third = Cmplx::new(Real::one() / Real::from_int(3), Real::zero());
                    x.ln().scale(&third).exp()
                }
            }

            // --- Special: Sinc = sin(x) / x ---
            Self::Sinc => {
                if x.coeffs[0].is_zero() {
                    // sinc(0) = 1 − x²/6 + x⁴/120 − …  (even series)
                    let mut coeffs = vec![Cmplx::zero(); deg + 1];
                    coeffs[0] = Cmplx::one();
                    // The recurrence a_n=-a_(n-2)/(n*(n+1)) avoids forming factorials.
                    let mut previous = Cmplx::one();
                    for n in (2..=deg).step_by(2) {
                        let denominator = Cmplx::new(
                            Real::from_int(
                                IntType::try_from(n).expect("Jet degree fits the integer backend"),
                            )
                            .to_float()
                                * Real::from_int(
                                    IntType::try_from(n + 1)
                                        .expect("Jet degree fits the integer backend"),
                                )
                                .to_float(),
                            Real::zero(),
                        );
                        previous = previous.neg().div(&denominator);
                        coeffs[n] = previous.clone();
                    }
                    // Local series at 0 composed with x itself (x₀ = 0, so
                    // dx = x): handles non-simple inputs like 2ε correctly.
                    Self::compose_local(&Jet::new(coeffs), x)
                } else {
                    x.sin().div(x)
                }
            }

            // --- Var / Const / algebraic ---
            Self::Var => x.clone(),
            Self::Const(c) => Jet::constant(c.clone(), deg),
            Self::Neg(a) => a.eval_jet(x).neg(),
            Self::Inv(a) => a.eval_jet(x).inv(),
            Self::Add(args) => {
                let mut acc = Jet::zero(deg);
                for a in args {
                    acc = acc.add(&a.eval_jet(x));
                }
                acc
            }
            Self::Mul(args) => {
                let mut acc = Jet::one(deg);
                for a in args {
                    acc = acc.mul(&a.eval_jet(x));
                }
                acc
            }
            Self::Compose(outer, inner) => {
                let inner_jet = inner.eval_jet(x);
                outer.eval_jet(&inner_jet)
            }

            // --- Functions with non-circular derivatives ---
            //
            // These functions' derivative ASTs are Jet-compatible (no reference
            // back to the function itself), so the derivative formula works.
            Self::Asin
            | Self::Acos
            | Self::Atan
            | Self::Acot
            | Self::Asec
            | Self::Acsc
            | Self::Asinh
            | Self::Acosh
            | Self::Atanh
            | Self::Acoth
            | Self::Asech
            | Self::Acsch => self.noncircular_eval_jet(x),
        }
    }

    /// Compose a local Taylor jet (coefficients at `x.coeffs[0]`) with the
    /// displacement `dx = x − x₀`, yielding the correct jet for general inputs.
    ///
    /// Local coefficients cannot be used directly for a general input series.
    /// Composition preserves the displacement encoded by its higher orders.
    /// Horner evaluation with truncation to `x`'s degree.
    fn compose_local(local: &Jet, x: &Jet) -> Jet {
        let deg = x.degree();
        if deg == 0 || (x.coeffs[1] == Cmplx::one() && x.coeffs.iter().skip(2).all(Cmplx::is_zero))
        {
            return local.clone();
        }

        let mut dx = x.coeffs.clone();
        if let Some(c0) = dx.first_mut() {
            *c0 = Cmplx::zero();
        }
        let mut acc = vec![Cmplx::zero(); deg + 1];
        for c in local.coeffs.iter().rev() {
            let mut next = vec![Cmplx::zero(); deg + 1];
            for m in 0..=deg {
                let mut sum = Cmplx::zero();
                for i in 0..=m {
                    if m - i <= deg {
                        sum += &acc[i].mul(&dx[m - i]);
                    }
                }
                next[m] = sum;
            }
            next[0] = next[0].add(c);
            acc = next;
        }
        Jet::new(acc)
    }

    /// Real-branch cube-root Jet: `y₀ = cbrt(x₀)` with coefficients from `y³ = x`.
    ///
    /// Matches scalar `Real::cbrt` on real inputs (unlike the ln/exp path,
    /// which takes the principal complex branch). From `y³ = x`:
    /// `3·y₀²·yₙ = xₙ − Σ_{i+j+k=n, i,j,k<n} yᵢyⱼyₖ`.
    fn cbrt_real_jet(x: &Jet) -> Jet {
        let deg = x.degree();
        if x.coeffs.iter().all(Cmplx::is_zero) {
            return Jet::zero(deg);
        }
        let y0 = Cmplx::new(x.coeffs[0].0.cbrt(), Real::zero());
        if y0.0.is_zero() {
            // x₀ = 0 with nonzero higher terms: branch point, derivatives
            // singular. Propagate NaN (callers decline non-finite terms).
            let nan = Cmplx::new(Real::nan(), Real::nan());
            return Jet::new(vec![nan; deg + 1]);
        }
        // denom = 3·y₀²
        let three = Cmplx::new(Real::one() + Real::one() + Real::one(), Real::zero());
        let denom = three.mul(&y0).mul(&y0);
        let mut y = vec![Cmplx::zero(); deg + 1];
        y[0] = y0;
        for n in 1..=deg {
            let mut sub = Cmplx::zero();
            for i in 0..n {
                // j <= n - i keeps k = n - i - j in range; k == n (i = j = 0)
                // is the 3·y₀²·yₙ term solved for, hence excluded.
                for j in 0..=n - i {
                    let k = n - i - j;
                    if k < n {
                        sub = sub.add(&y[i].mul(&y[j]).mul(&y[k]));
                    }
                }
            }
            y[n] = x.coeffs[n].sub(&sub).div(&denom);
        }
        Jet::new(y)
    }

    /// Fallback for functions that lack a direct Jet recurrence but whose
    /// derivative AST IS Jet-compatible (no circular reference back to self).
    /// Uses derivative relation: f'(x) = df(x), then recovers f's Taylor
    /// coefficients from df's and x's via the shift identity:
    ///   (n)·coeffs[n] = Σ_{i=0}^{n-1} `df_jet`[i] · (n-i) · x.coeffs[n-i]
    fn noncircular_eval_jet(&self, x: &Jet) -> Jet {
        let deg = x.degree();
        let mut coeffs = vec![Cmplx::zero(); deg + 1];
        coeffs[0] = self.eval(&x.coeffs[0]);

        if deg == 0 {
            return Jet::new(coeffs);
        }

        let df = self.derivative();
        // A non-finite derivative cannot provide finite higher coefficients.
        if matches!(&df, Self::Const(c) if !c.0.is_finite() || !c.1.is_finite()) {
            return Jet::new(coeffs);
        }

        let df_jet = df.eval_jet(x);
        for (n, c) in coeffs.iter_mut().enumerate().skip(1) {
            let mut sum = Cmplx::zero();
            for i in 0..n {
                let weight = Cmplx::new(
                    Real::from_int(
                        IntType::try_from(n - i).expect("Jet degree fits the integer backend"),
                    ),
                    Real::zero(),
                );
                sum += &df_jet.coeffs[i].mul(&weight).mul(&x.coeffs[n - i]);
            }
            *c = sum.div(&Cmplx::new(
                Real::from_int(IntType::try_from(n).expect("Jet degree fits the integer backend")),
                Real::zero(),
            ));
        }
        Jet::new(coeffs)
    }
}

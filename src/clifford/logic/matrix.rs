use core::{
    cmp::Ordering,
    ops::{AddAssign, Neg},
};

use alloc::vec::Vec;

use crate::{
    traits::Numeric,
    types::{Complex, RationalType, Real},
};

use super::{BladeAlgebra, CliffordNumber, GeneratorSet};
// Complex number (a + bi) using Real for both components

#[derive(Debug, Clone, PartialEq)]
pub struct Cmplx(pub Real, pub Real);

impl Cmplx {
    pub const fn new(re: Real, im: Real) -> Self {
        Self(re, im)
    }
    pub const fn zero() -> Self {
        Self(Real::zero(), Real::zero())
    }
    pub const fn one() -> Self {
        Self(Real::one(), Real::zero())
    }
    pub const fn i() -> Self {
        Self(Real::zero(), Real::one())
    }

    pub fn conj(&self) -> Self {
        Self(self.0.clone(), -&self.1)
    }
    pub fn abs_sq(&self) -> Real {
        &self.0 * &self.0 + &self.1 * &self.1
    }
    pub fn scale(&self, scalar: &Real) -> Self {
        Self(&self.0 * scalar, &self.1 * scalar)
    }

    pub fn add(&self, rhs: &Self) -> Self {
        Self(&self.0 + &rhs.0, &self.1 + &rhs.1)
    }
    pub fn sub(&self, rhs: &Self) -> Self {
        Self(&self.0 - &rhs.0, &self.1 - &rhs.1)
    }
    pub fn mul(&self, rhs: &Self) -> Self {
        Self(
            &self.0 * &rhs.0 - &self.1 * &rhs.1,
            &self.0 * &rhs.1 + &self.1 * &rhs.0,
        )
    }
    pub fn neg(&self) -> Self {
        Self(-&self.0, -&self.1)
    }

    fn typed(&self) -> Complex {
        Complex::from_parts(self.0.clone(), self.1.clone())
    }
    pub fn abs(&self) -> Real {
        self.typed().abs()
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_zero() && self.1.is_zero()
    }
    pub fn div(&self, rhs: &Self) -> Self {
        let value = self.typed().div(&rhs.typed());
        Self(value.re, value.im)
    }
    pub fn sqrt(&self) -> Self {
        let value = Numeric::sqrt(&self.typed());
        Self(value.re, value.im)
    }
    pub fn exp(&self) -> Self {
        let value = Numeric::exp(&self.typed());
        Self(value.re, value.im)
    }
    pub fn expm1(&self) -> Self {
        let value = Numeric::expm1(&self.typed());
        Self(value.re, value.im)
    }
    pub fn ln(&self) -> Self {
        let value = Numeric::ln(&self.typed());
        Self(value.re, value.im)
    }
    pub fn log1p(&self) -> Self {
        let value = Numeric::log1p(&self.typed());
        Self(value.re, value.im)
    }
    pub fn sin(&self) -> Self {
        let value = Numeric::sin(&self.typed());
        Self(value.re, value.im)
    }
    pub fn cos(&self) -> Self {
        let value = Numeric::cos(&self.typed());
        Self(value.re, value.im)
    }
    pub fn sinh(&self) -> Self {
        let value = Numeric::sinh(&self.typed());
        Self(value.re, value.im)
    }
    pub fn cosh(&self) -> Self {
        let value = Numeric::cosh(&self.typed());
        Self(value.re, value.im)
    }
    pub fn asin(&self) -> Self {
        let value = Numeric::asin(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acos(&self) -> Self {
        let value = Numeric::acos(&self.typed());
        Self(value.re, value.im)
    }
    pub fn atan(&self) -> Self {
        let value = Numeric::atan(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acot(&self) -> Self {
        let value = Numeric::acot(&self.typed());
        Self(value.re, value.im)
    }
    pub fn asec(&self) -> Self {
        let value = Numeric::asec(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acsc(&self) -> Self {
        let value = Numeric::acsc(&self.typed());
        Self(value.re, value.im)
    }
    pub fn asinh(&self) -> Self {
        let value = Numeric::asinh(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acosh(&self) -> Self {
        let value = Numeric::acosh(&self.typed());
        Self(value.re, value.im)
    }
    pub fn atanh(&self) -> Self {
        let value = Numeric::atanh(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acoth(&self) -> Self {
        let value = Numeric::acoth(&self.typed());
        Self(value.re, value.im)
    }
    pub fn asech(&self) -> Self {
        let value = Numeric::asech(&self.typed());
        Self(value.re, value.im)
    }
    pub fn acsch(&self) -> Self {
        let value = Numeric::acsch(&self.typed());
        Self(value.re, value.im)
    }
    pub fn cbrt(&self) -> Self {
        let value = Numeric::cbrt(&self.typed());
        Self(value.re, value.im)
    }
}
impl AddAssign<&Self> for Cmplx {
    fn add_assign(&mut self, rhs: &Self) {
        self.0 += &rhs.0;
        self.1 += &rhs.1;
    }
}
// 2×2 complex matrix [[a, b], [c, d]]

#[derive(Debug, Clone, PartialEq)]
pub struct Mat2C {
    pub a: Cmplx,
    pub b: Cmplx,
    pub c: Cmplx,
    pub d: Cmplx,
}

impl Mat2C {
    pub const fn new(a: Cmplx, b: Cmplx, c: Cmplx, d: Cmplx) -> Self {
        Self { a, b, c, d }
    }

    pub fn zero() -> Self {
        let z = Cmplx::zero();
        Self {
            a: z.clone(),
            b: z.clone(),
            c: z.clone(),
            d: z,
        }
    }

    pub const fn identity() -> Self {
        Self {
            a: Cmplx::one(),
            b: Cmplx::zero(),
            c: Cmplx::zero(),
            d: Cmplx::one(),
        }
    }

    pub fn det(&self) -> Cmplx {
        self.a.mul(&self.d).sub(&self.b.mul(&self.c))
    }
    pub fn trace(&self) -> Cmplx {
        self.a.add(&self.d)
    }

    /// Matrix inverse via the adjugate formula, or `None` when singular.
    ///
    /// Singularity is scale-aware (`|det| <= eps * (|ad| + |bc|)`), so
    /// exactly-singular integer matrices and numerically singular float
    /// matrices both decline. No eigenvectors involved.
    pub fn inverse(&self) -> Option<Self> {
        let ad = self.a.mul(&self.d);
        let bc = self.b.mul(&self.c);
        let det = ad.sub(&bc);
        let scale = ad.abs() + bc.abs();
        let tol = Real::epsilon() * scale;
        if det.abs().total_cmp(&tol) != Ordering::Greater {
            return None;
        }
        let inv = Cmplx::one().div(&det);
        Some(Self::new(
            self.d.mul(&inv),
            self.b.neg().mul(&inv),
            self.c.neg().mul(&inv),
            self.a.mul(&inv),
        ))
    }

    pub fn add(&self, rhs: &Self) -> Self {
        Self {
            a: self.a.add(&rhs.a),
            b: self.b.add(&rhs.b),
            c: self.c.add(&rhs.c),
            d: self.d.add(&rhs.d),
        }
    }

    pub fn mul(&self, rhs: &Self) -> Self {
        Self {
            a: self.a.mul(&rhs.a).add(&self.b.mul(&rhs.c)),
            b: self.a.mul(&rhs.b).add(&self.b.mul(&rhs.d)),
            c: self.c.mul(&rhs.a).add(&self.d.mul(&rhs.c)),
            d: self.c.mul(&rhs.b).add(&self.d.mul(&rhs.d)),
        }
    }

    pub fn scale(&self, scalar: &Cmplx) -> Self {
        Self {
            a: self.a.mul(scalar),
            b: self.b.mul(scalar),
            c: self.c.mul(scalar),
            d: self.d.mul(scalar),
        }
    }

    pub fn dagger(&self) -> Self {
        Self {
            a: self.a.conj(),
            b: self.c.conj(),
            c: self.b.conj(),
            d: self.d.conj(),
        }
    }

    /// Eigenvalues `(λ₁, λ₂)` via `(tr ± √(tr² − 4·det)) / 2`.
    pub fn eigenvalues(&self) -> (Cmplx, Cmplx) {
        let tr = self.trace();
        let det = self.det();
        let four = Cmplx::new(
            Real::one() + Real::one() + Real::one() + Real::one(),
            Real::zero(),
        );
        let disc = tr.mul(&tr).sub(&det.mul(&four));
        let sqrt_disc = disc.sqrt();
        let half = Cmplx::new(Real::one() / (Real::one() + Real::one()), Real::zero());
        (tr.add(&sqrt_disc).mul(&half), tr.sub(&sqrt_disc).mul(&half))
    }
}

/// Small-algebra matrix embedding: generator-set embeddability checks,
/// blade-matrix construction, multivector projection, and 2x2
/// eigendecomposition.
pub struct SmallEmbed;

impl SmallEmbed {
    /// Maximum generators that can be embedded in Mat(2,ℂ).
    const MAX_EMBED_GENS: usize = 3;

    /// Check whether `gens` can be embedded in Mat(2,ℂ) (requires `gens.len() + num_nilpotents <= 3`).
    ///
    /// This is a necessary budget check, not a fidelity guarantee: some
    /// signatures within budget still collapse blades (e.g. `[1, 1, -1]`
    /// maps `e1` and `-e23` to the same image). The definitive round-trip
    /// check lives in [`Self::compute_blade_basis`], which returns `None`
    /// for unfaithful signatures so callers fall back to the large path.
    pub fn can_embed(gens: &GeneratorSet) -> bool {
        let mut num_nilpotent = 0;
        for i in 0..gens.len() {
            if gens.metric_at(i) == 0 {
                num_nilpotent += 1;
            }
        }
        gens.len() + num_nilpotent <= Self::MAX_EMBED_GENS
    }

    /// Return the 2×2 Pauli matrix for a given generator position.
    /// Position 0 -> `sigma_x`, 1 -> `sigma_y`, 2 -> `sigma_z`.
    fn pauli_matrix(pos: usize) -> Mat2C {
        debug_assert!(
            pos < 3,
            "Pauli matrix position must be 0, 1, or 2, got {pos}"
        );
        let z = Cmplx::zero();
        let one = Cmplx::one();
        let i = Cmplx::i();
        if pos == 0 {
            Mat2C::new(z.clone(), one.clone(), one, z)
        } else if pos == 1 {
            Mat2C::new(z.clone(), i.neg(), i, z)
        } else {
            Mat2C::new(one.clone(), z.clone(), z, one.neg())
        }
    }

    /// Compute ALL blade matrices for a given generator set.
    /// Returns `None` if the generator set is too large.
    pub fn compute_blade_basis(gens: &GeneratorSet) -> Option<Vec<Mat2C>> {
        if !Self::can_embed(gens) {
            return None;
        }
        let n = gens.len();
        let mut gen_mats: Vec<Mat2C> = Vec::with_capacity(n);
        let mut pauli_idx = 0;
        for pos in 0..n {
            let metric = gens.metric_at(pos);
            if metric == 1 {
                gen_mats.push(Self::pauli_matrix(pauli_idx));
                pauli_idx += 1;
            } else if metric == -1 {
                gen_mats.push(Self::pauli_matrix(pauli_idx).scale(&Cmplx::i()));
                pauli_idx += 1;
            } else {
                let g1 = Self::pauli_matrix(pauli_idx);
                let g2 = Self::pauli_matrix(pauli_idx + 1);
                let i_g2 = g2.scale(&Cmplx::i());
                let mut n_mat = g1.add(&i_g2);
                let half = Cmplx::new(
                    Real::from_rational(RationalType::from_parts(1, 2)),
                    Real::zero(),
                );
                n_mat = n_mat.scale(&half);
                gen_mats.push(n_mat);
                pauli_idx += 2;
            }
        }

        if gens.complex_blade().is_some() {
            let mut product = Mat2C::identity();
            for generator in &gen_mats[..n - 1] {
                product = product.mul(generator);
            }
            let square = product.mul(&product).a;
            gen_mats[n - 1] = product.scale(&Cmplx::i().div(&square));
        }
        let count = 1_usize << n;
        let mut basis = Vec::with_capacity(count);
        for blade in 0..count {
            let mut m = Mat2C::identity();
            for (i, gm) in gen_mats.iter().enumerate() {
                if (blade >> i) & 1 == 1 {
                    m = m.mul(gm);
                }
            }
            basis.push(m);
        }
        // Fidelity gate: the trace projection must round-trip every blade
        // exactly. Some in-budget signatures collapse distinct blades onto
        // the same image (e.g. `[1, 1, -1]`); those return `None` so callers
        // fall back to the faithful large path instead of silently leaking
        // coefficients across blades.
        for (blade, bm) in basis.iter().enumerate() {
            let coeffs = Self::project_onto_basis(bm, &basis);
            for (b, c) in coeffs.iter().enumerate() {
                let ok = if b == blade { c.is_one() } else { c.is_zero() };
                if !ok {
                    return None;
                }
            }
        }
        Some(basis)
    }

    /// Convert a [`CliffordNumber`] to a 2×2 complex matrix.
    /// Returns `None` if the generator set cannot be embedded in Mat(2,ℂ).
    pub fn to_matrix_with_basis(mv: &CliffordNumber, blade_mats: &[Mat2C]) -> Mat2C {
        let mut result = Mat2C::zero();
        for (blade, bm) in blade_mats.iter().enumerate() {
            let coeff = mv.coeff(blade);
            if coeff.is_zero() {
                continue;
            }
            let scaled = bm.scale(&Cmplx::new(coeff.clone(), Real::zero()));
            result = result.add(&scaled);
        }
        result
    }

    /// Project a Mat(2,ℂ) matrix back onto the blade basis of `gens`, returning
    /// the blade coefficients.
    fn project_onto_basis(mat: &Mat2C, blade_mats: &[Mat2C]) -> Vec<Real> {
        blade_mats
            .iter()
            .map(|bm| {
                let prod = mat.mul(&bm.dagger());
                let numer = prod.trace().0;
                let norm_bm = bm.mul(&bm.dagger());
                let denom = norm_bm.trace().0;
                if denom.is_zero() {
                    Real::zero()
                } else {
                    numer / denom
                }
            })
            .collect()
    }

    /// Convert a 2×2 complex matrix back to a [`CliffordNumber`] for a given
    /// generator set. The generator set must have been used to produce the matrix
    /// (same blade basis).
    pub fn from_matrix_with_gens(
        mat: &Mat2C,
        gens: &GeneratorSet,
        blade_mats: &[Mat2C],
    ) -> CliffordNumber {
        let coeffs = Self::project_onto_basis(mat, blade_mats);

        let count = BladeAlgebra::coeff_count_unchecked(
            u8::try_from(gens.len()).expect("can_embed ensures n <= 3 fits in u8"),
        );
        let mut mv = CliffordNumber::zero_unchecked(gens.clone());
        for (blade, c) in coeffs.iter().enumerate().take(count) {
            if !c.is_zero() {
                mv.coeffs_mut_slice()[blade] = c.clone();
            }
        }
        mv
    }
    // Eigendecomposition

    /// Eigenvector for `lambda`: solves `(M - λI)v = 0`.
    ///
    /// Selects the larger-norm candidate from (b, lambda-a) and (lambda-d, c).
    /// Returns None when both vectors vanish.
    fn eigenvector_for(mat: &Mat2C, lambda: &Cmplx) -> Option<(Cmplx, Cmplx)> {
        let row1 = (mat.b.clone(), lambda.sub(&mat.a));
        let row2 = (lambda.sub(&mat.d), mat.c.clone());
        let n1 = row1.0.abs_sq() + row1.1.abs_sq();
        let n2 = row2.0.abs_sq() + row2.1.abs_sq();
        if n1.is_zero() && n2.is_zero() {
            return None;
        }
        if n1.total_cmp(&n2) == Ordering::Less {
            Some(row2)
        } else {
            Some(row1)
        }
    }

    /// Returns (U, U^-1, lambda1, lambda2) for a diagonalizable 2x2 matrix.
    /// Returns None for a repeated non-diagonal eigenvalue or singular U.
    pub fn eigendecompose(mat: &Mat2C) -> Option<(Mat2C, Mat2C, Cmplx, Cmplx)> {
        // For diagonal input, U = I and eigenvalues retain their diagonal order.
        if mat.b.is_zero() && mat.c.is_zero() {
            return Some((
                Mat2C::identity(),
                Mat2C::identity(),
                mat.a.clone(),
                mat.d.clone(),
            ));
        }

        let (lambda1, lambda2) = mat.eigenvalues();

        // Check for defective case: λ₁ ≈ λ₂ but A ≠ λI
        let is_defective = lambda1.sub(&lambda2).abs_sq().is_zero();
        if is_defective {
            return None; // not diagonalizable
        }

        // Eigenvector for λ₁: solve (M - λ₁I)v = 0
        // Using the row with larger norm for numerical stability
        let (v1_a, v1_b) = Self::eigenvector_for(mat, &lambda1)?;

        // Eigenvector for λ₂: same approach
        let (v2_a, v2_b) = Self::eigenvector_for(mat, &lambda2)?;

        // U = [v1 v2] (columns are eigenvectors)
        let u = Mat2C::new(v1_a, v2_a, v1_b, v2_b);

        // U⁻¹ = (1/det) · [[d, -b], [-c, a]]
        let det_u = u.det();
        // For nonzero det(U), its reciprocal is conj(det(U)) / |det(U)|^2.
        let det_abs_sq = det_u.abs_sq();
        // A singular eigenvector matrix cannot define the decomposition.
        if det_abs_sq.is_zero() {
            return None;
        }
        let inv_det = Cmplx::new(det_u.0.clone(), det_u.1.neg()).scale(&(Real::one() / det_abs_sq));

        let u_inv = Mat2C::new(
            u.d.mul(&inv_det),
            u.b.neg().mul(&inv_det),
            u.c.neg().mul(&inv_det),
            u.a.mul(&inv_det),
        );

        Some((u, u_inv, lambda1, lambda2))
    }
}

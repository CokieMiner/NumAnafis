use core::{cmp::Ordering, ops::Add};

use crate::{
    traits::Numeric,
    types::{Complex, RationalType, Real},
};

use super::{
    CliffordNumber, Cmplx, EigenSolver, Jet, LargeEmbed, LargeSpectral, Mat2C, SmallEmbed,
    SpectralFn,
};
// SpectralDispatch — backend selection for spectral evaluation

/// Scalar complex evaluation, faithful 2x2 embedding, and general matrix evaluation.
pub struct SpectralDispatch;

impl SpectralDispatch {
    pub fn apply_via_spectral(mv: &CliffordNumber, ast: &SpectralFn) -> CliffordNumber {
        let gens = mv.generator_set();

        if let Some((re, im)) = mv.central_scalar_parts() {
            let value = ast.eval(&Cmplx::new(re.clone(), im.clone()));
            return CliffordNumber::scalar(gens.clone(), Complex::from_parts(value.0, value.1))
                .unwrap_or_else(|_| CliffordNumber::nan(gens.clone()));
        }

        // Path 1: ≤3 generators via Mat(2,ℂ) spectral decomposition
        if let Some(blade_mats) = SmallEmbed::compute_blade_basis(gens) {
            let mat = SmallEmbed::to_matrix_with_basis(mv, &blade_mats);

            // Standard Eigendecomposition
            if let Some((u, u_inv, l1, l2)) = SmallEmbed::eigendecompose(&mat) {
                let fd = Mat2C::new(ast.eval(&l1), Cmplx::zero(), Cmplx::zero(), ast.eval(&l2));
                let result_mat = u.mul(&fd).mul(&u_inv);
                let projected = SmallEmbed::from_matrix_with_gens(&result_mat, gens, &blade_mats);
                if Self::projection_valid(&result_mat, &projected, &blade_mats) {
                    return projected;
                }
                return CliffordNumber::nan(gens.clone());
            }

            // Defective matrix fallback: compute via Schur form T = Q^H * A * Q
            let (q_mat, t_mat) = Self::schur_decompose_2x2(&mat);
            let lambda = t_mat.a.clone();
            let f_lambda = ast.eval(&lambda);
            // First-order jets evaluate the derivative and removable singularities
            // such as sinc'(0) = 0.
            let jet1 = ast.eval_jet(&Jet::variable(&lambda, 1));
            let Some(f_prime_lambda) = jet1.coeffs.get(1).cloned() else {
                return CliffordNumber::nan(gens.clone());
            };
            if f_lambda.0.is_finite() && f_prime_lambda.0.is_finite() {
                let fd = Mat2C::new(
                    f_lambda.clone(),
                    t_mat.b.mul(&f_prime_lambda),
                    Cmplx::zero(),
                    f_lambda,
                );
                let q_dagger = q_mat.dagger();
                let result_mat = q_mat.mul(&fd).mul(&q_dagger);
                let projected = SmallEmbed::from_matrix_with_gens(&result_mat, gens, &blade_mats);
                if Self::projection_valid(&result_mat, &projected, &blade_mats) {
                    return projected;
                }
                return CliffordNumber::nan(gens.clone());
            }
            return CliffordNumber::nan(gens.clone());
        }

        // Path 2: Fallback to large_mat.rs for large or highly nilpotent algebras
        LargeSpectral::large_apply_via_spectral(mv, ast)
    }

    pub fn extract_scalar(mv: &CliffordNumber) -> Real {
        let blade0 = mv.coeff(0);
        if blade0.is_zero() {
            Real::zero()
        } else {
            blade0.clone()
        }
    }

    pub fn is_pure_scalar(mv: &CliffordNumber) -> bool {
        for blade in 1..mv.blade_count() {
            if !mv.coeff(blade).is_zero() {
                return false;
            }
        }
        true
    }

    // Closure-based functional calculus for binary evaluator operations.
    pub fn apply_via_spectral_closure<F>(mv: &CliffordNumber, f: &F) -> CliffordNumber
    where
        F: Fn(&Cmplx) -> Cmplx,
    {
        let gens = mv.generator_set();
        if let Some((re, im)) = mv.central_scalar_parts() {
            let value = f(&Cmplx::new(re.clone(), im.clone()));
            return CliffordNumber::scalar(gens.clone(), Complex::from_parts(value.0, value.1))
                .unwrap_or_else(|_| CliffordNumber::nan(gens.clone()));
        }

        if let Some(blade_mats) = SmallEmbed::compute_blade_basis(gens) {
            let mat = SmallEmbed::to_matrix_with_basis(mv, &blade_mats);
            if let Some((u, u_inv, l1, l2)) = SmallEmbed::eigendecompose(&mat) {
                let fd = Mat2C::new(f(&l1), Cmplx::zero(), Cmplx::zero(), f(&l2));
                let result_mat = u.mul(&fd).mul(&u_inv);
                let projected = SmallEmbed::from_matrix_with_gens(&result_mat, gens, &blade_mats);
                if Self::projection_valid(&result_mat, &projected, &blade_mats) {
                    return projected;
                }
                return CliffordNumber::nan(gens.clone());
            }
            // Defective matrix fallback via Schur form + numeric derivative
            let (q_mat, t_mat) = Self::schur_decompose_2x2(&mat);
            let lambda = t_mat.a.clone();
            let f_lambda = f(&lambda);
            if f_lambda.0.is_finite() && f_lambda.1.is_zero() {
                // Numeric derivative via central finite differences with a
                // scale-aware step h ≈ eps^(1/3)·max(1,|λ|). A bare-epsilon
                // step rounds away below float resolution for |λ| ≳ 1,
                // silently zeroing the derivative.
                let lam_abs = lambda.abs();
                let one = Real::one();
                let scale = if lam_abs.total_cmp(&one) == Ordering::Greater {
                    lam_abs
                } else {
                    one
                };
                let h = Real::epsilon().cbrt() * scale;
                let h_cmplx = Cmplx::new(h.clone(), Real::zero());
                let f_plus = f(&lambda.add(&h_cmplx));
                let f_minus = f(&lambda.sub(&h_cmplx));
                let two_h = Cmplx::new(Real::from_int(2) * h, Real::zero());
                let f_prime_lambda = (f_plus.sub(&f_minus)).div(&two_h);
                if f_prime_lambda.0.is_finite() && f_prime_lambda.1.is_zero() {
                    let fd = Mat2C::new(
                        f_lambda.clone(),
                        t_mat.b.mul(&f_prime_lambda),
                        Cmplx::zero(),
                        f_lambda,
                    );
                    let q_dagger = q_mat.dagger();
                    let result_mat = q_mat.mul(&fd).mul(&q_dagger);
                    let projected =
                        SmallEmbed::from_matrix_with_gens(&result_mat, gens, &blade_mats);
                    if Self::projection_valid(&result_mat, &projected, &blade_mats) {
                        return projected;
                    }
                    return CliffordNumber::nan(gens.clone());
                }
            }
            return CliffordNumber::nan(gens.clone());
        }
        // Fallback to large_mat.rs
        LargeSpectral::large_apply_via_spectral_closure(mv, f)
    }

    /// Verify a projected multivector re-embeds to its source matrix.
    ///
    /// Spectral images need not lie in the represented algebra (e.g. √e ∉
    /// Cl(1,0)); the trace projection then returns a real multivector that
    /// does not satisfy the function equation. Callers decline (NaN) when
    /// the residual exceeds tolerance instead of returning it.
    fn projection_valid(mat: &Mat2C, mv: &CliffordNumber, blade_mats: &[Mat2C]) -> bool {
        let rt = SmallEmbed::to_matrix_with_basis(mv, blade_mats);
        let mut scale = Real::zero();
        for entry in [&mat.a, &mat.b, &mat.c, &mat.d] {
            let a = entry.abs();
            if a.total_cmp(&scale) == Ordering::Greater {
                scale = a;
            }
        }
        let tol = Real::epsilon() * scale * Real::from_int(100);
        for (x, y) in [&mat.a, &mat.b, &mat.c, &mat.d]
            .into_iter()
            .zip([&rt.a, &rt.b, &rt.c, &rt.d])
        {
            if x.sub(y).abs().total_cmp(&tol) == Ordering::Greater {
                return false;
            }
        }
        true
    }

    /// Compute the Schur decomposition of a 2×2 matrix.
    /// Returns `(Q, T)` where `Q` is unitary and `T = Q^H * A * Q` is upper triangular.
    /// For a defective matrix, T = [[λ, x], [0, λ]] with repeated eigenvalue λ.
    fn schur_decompose_2x2(mat: &Mat2C) -> (Mat2C, Mat2C) {
        // Eigenvalue: λ = trace / 2 (for 2x2, repeated eigenvalue case)
        let trace = mat.a.add(&mat.d);
        let half = Cmplx::new(
            Real::from_rational(RationalType::from_parts(1, 2)),
            Real::zero(),
        );
        let lambda = trace.mul(&half);

        // Use the larger of [b, lambda-a] and [lambda-d, c].
        // A lower Jordan block has eigenvector [0,1]; fixing v1=1 fails.
        let first = (mat.b.clone(), lambda.sub(&mat.a));
        let second = (lambda.sub(&mat.d), mat.c.clone());
        let first_norm = first.0.abs_sq() + first.1.abs_sq();
        let second_norm = second.0.abs_sq() + second.1.abs_sq();
        let (v1, v2) = if first_norm.is_zero() && second_norm.is_zero() {
            (Cmplx::one(), Cmplx::zero())
        } else if first_norm >= second_norm {
            first
        } else {
            second
        };

        // Normalize v
        let v_norm = v1.abs_sq().add(&v2.abs_sq()).sqrt();
        let v_norm_c = Cmplx::new(v_norm, Real::zero());
        let v1n = v1.div(&v_norm_c);
        let v2n = v2.div(&v_norm_c);

        // Find orthogonal vector w = [-conj(v2), conj(v1)]
        let w1 = v2n.conj().neg();
        let w2 = v1n.conj();

        // Build Q = [v, w]
        let q_mat = Mat2C::new(v1n, w1, v2n, w2);

        // Compute T = Q^H * A * Q
        let qh = q_mat.dagger();
        let qh_a = qh.mul(mat);
        let t_mat = qh_a.mul(&q_mat);

        (q_mat, t_mat)
    }
}
// Geometric inverse via spectral decomposition

impl CliffordNumber {
    /// Computes the geometric inverse of the Clifford number.
    #[must_use]
    pub fn geometric_inverse(&self) -> Self {
        let gens = self.generator_set();

        // Path 1: ≤3 generators via Mat(2,ℂ) adjugate inverse.
        // Adjugate inversion with a scale-aware singularity test.
        if let Some(blade_mats) = SmallEmbed::compute_blade_basis(gens) {
            let mat = SmallEmbed::to_matrix_with_basis(self, &blade_mats);
            if let Some(inv) = mat.inverse() {
                return SmallEmbed::from_matrix_with_gens(&inv, gens, &blade_mats);
            }
        }

        // Path 2: LU solve on the large-matrix embedding. Also
        // eigenvector-free; declines (numerically) singular inputs.
        if let Some(basis) = LargeEmbed::blade_basis(gens) {
            let mat = LargeEmbed::to_matrix(self, &basis);
            if let Some(inv) = EigenSolver::lu_solve(&mat) {
                return LargeEmbed::from_matrix(&inv, gens, &basis);
            }
        }
        Self::nan(gens.clone())
    }
}

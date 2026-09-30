//! General-dimension eigensolver for complex matrices.
//!
//! Hessenberg reduction with QR iteration, Schur decomposition, and
//! eigenvector recovery from the Schur form. Backs spectral evaluation
//! for algebras beyond the 2x2 embedding range.

use core::cmp::Ordering;

use alloc::{vec, vec::Vec};

use crate::{traits::Numeric, types::Real};

use super::{Cmplx, MatC};
// Heuristic native work cap; no general QR convergence bound is asserted.
// TODO(mp): scale this work limit with precision and convergence diagnostics.
// Iteration exhaustion must continue to decline the decomposition.
const QR_ITERATIONS_PER_DIMENSION: usize = 200;

// Eigendecomposition (QR Algorithm)

fn givens_rotation(x_val: &Cmplx, y_val: &Cmplx) -> (Cmplx, Cmplx) {
    if y_val.abs_sq().is_zero() {
        return (Cmplx::one(), Cmplx::zero());
    }
    if x_val.abs_sq().is_zero() {
        let mag = y_val.abs_sq().sqrt();
        let inv_mag = Cmplx::one().div(&Cmplx::new(mag, Real::zero()));
        return (Cmplx::zero(), y_val.conj().mul(&inv_mag));
    }
    let r_val = (x_val.abs_sq() + y_val.abs_sq()).sqrt();
    let x_mag = x_val.abs_sq().sqrt();

    let inv_r = Cmplx::one().div(&Cmplx::new(r_val, Real::zero()));
    let cos_val = Cmplx::new(x_mag.clone(), Real::zero()).mul(&inv_r);

    let inv_x_mag = Cmplx::one().div(&Cmplx::new(x_mag, Real::zero()));
    let sign_x = x_val.mul(&inv_x_mag);
    let sin_val = y_val.conj().mul(&sign_x).mul(&inv_r);

    (cos_val, sin_val)
}

fn hessenberg_reduce(a: &mut MatC, q: &mut MatC, buf_v: &mut Vec<Cmplx>) {
    let n = a.dim;
    if n <= 2 {
        return;
    }
    for k in 0..n - 2 {
        let mut norm_sq = Real::zero();
        for i in k + 1..n {
            norm_sq = norm_sq + a.get(i, k).abs_sq();
        }
        let norm = norm_sq.sqrt();
        if norm.is_zero() {
            continue;
        }

        buf_v.clear();
        for i in k + 1..n {
            buf_v.push(a.get(i, k).clone());
        }

        let x0 = &buf_v[0];
        let x0_mag = x0.abs_sq().sqrt();
        let phase = if x0_mag.is_zero() {
            Cmplx::one()
        } else {
            let inv_mag = Cmplx::one().div(&Cmplx::new(x0_mag, Real::zero()));
            x0.mul(&inv_mag)
        };
        let norm_c = Cmplx::new(norm, Real::zero());
        let shifted = phase.mul(&norm_c);
        buf_v[0] = buf_v[0].add(&shifted);

        let mut v_norm_sq = Cmplx::zero();
        for vi in buf_v.iter() {
            v_norm_sq = v_norm_sq.add(&vi.conj().mul(vi));
        }
        if v_norm_sq.abs_sq().is_zero() {
            continue;
        }

        let two = Cmplx::new(Real::one() + Real::one(), Real::zero());
        let two_over_v_norm_sq = two.div(&v_norm_sq);

        for j in k..n {
            let mut v_star_a = Cmplx::zero();
            for i in k + 1..n {
                v_star_a = v_star_a.add(&buf_v[i - k - 1].conj().mul(a.get(i, j)));
            }
            let scaled = v_star_a.mul(&two_over_v_norm_sq);
            for i in k + 1..n {
                let term = buf_v[i - k - 1].mul(&scaled);
                a.set(i, j, a.get(i, j).sub(&term));
            }
        }

        for i in k + 2..n {
            a.set(i, k, Cmplx::zero());
        }

        for i in 0..n {
            let mut a_v = Cmplx::zero();
            for j in k + 1..n {
                a_v = a_v.add(&a.get(i, j).mul(&buf_v[j - k - 1]));
            }
            let scaled = a_v.mul(&two_over_v_norm_sq);
            for j in k + 1..n {
                let term = scaled.mul(&buf_v[j - k - 1].conj());
                a.set(i, j, a.get(i, j).sub(&term));
            }
        }

        for i in 0..n {
            let mut q_v = Cmplx::zero();
            for j in k + 1..n {
                q_v = q_v.add(&q.get(i, j).mul(&buf_v[j - k - 1]));
            }
            let scaled = q_v.mul(&two_over_v_norm_sq);
            for j in k + 1..n {
                let term = scaled.mul(&buf_v[j - k - 1].conj());
                q.set(i, j, q.get(i, j).sub(&term));
            }
        }
    }
}

fn qr_step_hessenberg(
    a: &mut MatC,
    q_total: &mut MatC,
    active: usize,
    buf_givens: &mut Vec<(Cmplx, Cmplx)>,
) {
    let dim = a.dim;
    let mut shift = Cmplx::zero();
    if active >= 2 {
        let n1 = active - 1;
        let n2 = active - 2;
        let a_11 = a.get(n1, n1);
        let a_22 = a.get(n2, n2);
        let a_12 = a.get(n1, n2);
        let a_21 = a.get(n2, n1);

        let tr = a_22.add(a_11);
        let det = a_22.mul(a_11).sub(&a_21.mul(a_12));

        let four = Cmplx::new(
            Real::one() + Real::one() + Real::one() + Real::one(),
            Real::zero(),
        );
        let disc = tr.mul(&tr).sub(&four.mul(&det));
        let sqrt_disc = disc.sqrt();

        let inv_two = Cmplx::one().div(&Cmplx::new(Real::one() + Real::one(), Real::zero()));
        let l1 = tr.add(&sqrt_disc).mul(&inv_two);
        let l2 = tr.sub(&sqrt_disc).mul(&inv_two);

        if l1.sub(a_11).abs_sq().total_cmp(&l2.sub(a_11).abs_sq()) == Ordering::Less {
            shift = l1;
        } else {
            shift = l2;
        }
    }

    for i in 0..active {
        let val = a.get(i, i).sub(&shift);
        a.set(i, i, val);
    }

    buf_givens.clear();
    for i in 0..active - 1 {
        let x_val = a.get(i, i).clone();
        let y_val = a.get(i + 1, i).clone();
        let (cos_val, sin_val) = givens_rotation(&x_val, &y_val);
        buf_givens.push((cos_val.clone(), sin_val.clone()));

        for j in i..active {
            let a_ij = a.get(i, j).clone();
            let a_ip1j = a.get(i + 1, j).clone();
            a.set(i, j, cos_val.mul(&a_ij).add(&sin_val.mul(&a_ip1j)));
            a.set(
                i + 1,
                j,
                sin_val.conj().neg().mul(&a_ij).add(&cos_val.mul(&a_ip1j)),
            );
        }
    }

    for (i, (cos_val, sin_val)) in buf_givens.iter().cloned().enumerate() {
        for j in 0..=i + 1 {
            let a_ji = a.get(j, i).clone();
            let a_jip1 = a.get(j, i + 1).clone();
            a.set(j, i, a_ji.mul(&cos_val).add(&a_jip1.mul(&sin_val.conj())));
            a.set(
                j,
                i + 1,
                a_ji.mul(&sin_val.neg()).add(&a_jip1.mul(&cos_val)),
            );
        }
    }

    for i in 0..active {
        let val = a.get(i, i).add(&shift);
        a.set(i, i, val);
    }

    for (i, (cos_val, sin_val)) in buf_givens.iter().cloned().enumerate() {
        for j in 0..dim {
            let q_ji = q_total.get(j, i).clone();
            let q_jip1 = q_total.get(j, i + 1).clone();
            // Q accumulates the QR factor Q_fac = G^H (since R = G(A-σI)):
            // column i takes +scalar̄, column i+1 takes -scalar.
            q_total.set(j, i, q_ji.mul(&cos_val).add(&q_jip1.mul(&sin_val.conj())));
            q_total.set(
                j,
                i + 1,
                q_ji.mul(&sin_val.neg()).add(&q_jip1.mul(&cos_val)),
            );
        }
    }
}
// Eigenvector recovery from Schur form.

/// Check if an upper triangular matrix is numerically diagonal.
fn is_triangular_diagonal(t: &MatC, tol: &Real) -> bool {
    for i in 0..t.dim {
        for j in i + 1..t.dim {
            if t.get(i, j).abs().total_cmp(tol) == Ordering::Greater {
                return false;
            }
        }
    }
    true
}

/// Check whether all eigenvalues are numerically distinct.
fn eigenvalues_are_distinct(eigvals: &[Cmplx], tol: &Real) -> bool {
    for i in 0..eigvals.len() {
        for j in i + 1..eigvals.len() {
            if eigvals[i].sub(&eigvals[j]).abs().total_cmp(tol) == Ordering::Less {
                return false;
            }
        }
    }
    true
}

/// Compute eigenvectors of an upper triangular matrix with **distinct** eigenvalues.
///
/// Returns `V` where each column is an eigenvector, normalized so that `V[i][i] = 1`.
/// `V` is upper triangular with unit diagonal.
fn eigenvectors_from_triangular_distinct(t: &MatC) -> Option<MatC> {
    let n = t.dim;
    let eigvals = EigenSolver::eigenvalues_from_triangular(t);
    let tol = Real::epsilon();

    if !eigenvalues_are_distinct(&eigvals, &tol) {
        return None;
    }

    let mut v = MatC::identity(n); // V[:, i] = eigenvector for λ_i

    for (i, eig_i) in eigvals.iter().enumerate() {
        // For eigenvector v_i (column i), components > i are 0, component i is 1.
        // Back-substitute for components j = i-1, i-2, ..., 0:
        //   v_i[j] = -sum_{k=j+1}^{i} T[j][k] * v_i[k] / (T[j][j] - λ_i)
        for j in (0..i).rev() {
            let denom = t.get(j, j).sub(eig_i);
            let mut sum = Cmplx::zero();
            for k in j + 1..=i {
                sum += &t.get(j, k).mul(v.get(k, i));
            }
            v.set(j, i, sum.neg().div(&denom));
        }
    }
    Some(v)
}

/// Compute the inverse of an upper triangular matrix with unit diagonal.
/// Result is also upper triangular with unit diagonal.
fn inverse_upper_triangular(v: &MatC) -> MatC {
    let n = v.dim;
    let mut vinv = MatC::identity(n);

    // For i < j: (V^{-1})[i][j] = -sum_{k=i+1}^{j} V[i][k] * (V^{-1})[k][j]
    for j in 1..n {
        for i in (0..j).rev() {
            let mut sum = Cmplx::zero();
            for k in i + 1..=j {
                sum += &v.get(i, k).mul(vinv.get(k, j));
            }
            vinv.set(i, j, sum.neg());
        }
    }
    vinv
}

/// Numerical rank from the spectrum of the Gram matrix B^H B.
///
/// Counts Gram eigenvalues above tol times the largest Gram eigenvalue.
/// Equivalently, singular values exceed sqrt(tol) times the largest singular value.
/// Schur nonconvergence returns zero so callers reject the attempted decomposition.
fn rank_upper_triangular(mat: &MatC, tol: &Real) -> usize {
    let gram = mat.dagger().mul(mat);
    let Some((_, t)) = EigenSolver::schur_decomposition(&gram) else {
        return 0;
    };
    let mut smax = Real::zero();
    for i in 0..t.dim {
        let scalar = t.get(i, i).abs();
        if scalar.total_cmp(&smax) == Ordering::Greater {
            smax = scalar;
        }
    }
    let cutoff = tol * &smax;
    let mut scalar_ratio = 0;
    for i in 0..t.dim {
        if t.get(i, i).abs().total_cmp(&cutoff) == Ordering::Greater {
            scalar_ratio += 1;
        }
    }
    scalar_ratio
}
// EigenSolver — public eigensolver entry points

/// General-dimension eigensolver: Schur decomposition, eigenvalue
/// extraction, and eigenvector recovery from the Schur form.
pub struct EigenSolver;

impl EigenSolver {
    pub fn schur_decomposition(mat: &MatC) -> Option<(MatC, MatC)> {
        let mut a = mat.clone();
        let mut q = MatC::identity(mat.dim);

        // Pre-allocate buffers to avoid repeated heap allocations in the QR loop
        let mut buf_v: Vec<Cmplx> = Vec::with_capacity(mat.dim);
        let mut buf_givens: Vec<(Cmplx, Cmplx)> = Vec::with_capacity(mat.dim);

        let mut active = mat.dim;
        // Relative deflation tolerance, scaled by neighboring diagonal magnitudes.
        let tol = Real::epsilon();

        hessenberg_reduce(&mut a, &mut q, &mut buf_v);

        for _ in 0..mat.dim.checked_mul(QR_ITERATIONS_PER_DIMENSION)? {
            if active <= 1 {
                break;
            }

            let mut deflated = false;
            let sub_diag_abs = a.get(active - 1, active - 2).abs();

            let diag1 = a.get(active - 1, active - 1).abs();
            let diag2 = a.get(active - 2, active - 2).abs();
            let scale = diag1 + diag2;
            // Purely relative tolerance: no absolute floor, so tiny matrices
            // are judged against their own scale. An exactly-zero subdiagonal
            // always deflates; anything else must beat eps*100 relatively.
            let hundred = Real::from_int(100);
            let effective_tol = &(&tol * &scale) * &hundred;

            if sub_diag_abs.is_zero() || sub_diag_abs.total_cmp(&effective_tol) == Ordering::Less {
                a.set(active - 1, active - 2, Cmplx::zero());
                active -= 1;
                deflated = true;
            }

            if !deflated {
                qr_step_hessenberg(&mut a, &mut q, active, &mut buf_givens);
            }
        }

        if active > 1 {
            return None;
        }

        Some((q, a))
    }

    /// Extract eigenvalues from the diagonal of an upper triangular matrix.
    pub fn eigenvalues_from_triangular(t: &MatC) -> Vec<Cmplx> {
        let n = t.dim;
        let mut eigvals = Vec::with_capacity(n);
        for i in 0..n {
            eigvals.push(t.get(i, i).clone());
        }
        eigvals
    }

    /// Attempt a full eigendecomposition from a Schur factor `T`.
    ///
    /// Returns `(V, V_inv)` where `T = V * D * V^{-1}`,
    /// or `None` if the matrix is defective.
    pub fn eigendecompose_from_schur(t: &MatC) -> Option<(MatC, MatC)> {
        let tol = Real::epsilon();

        // Fast path: T is already diagonal → V = I
        if is_triangular_diagonal(t, &tol) {
            return Some((MatC::identity(t.dim), MatC::identity(t.dim)));
        }

        // Try to compute eigenvectors assuming distinct eigenvalues (most common case)
        if let Some(v) = eigenvectors_from_triangular_distinct(t) {
            let vinv = inverse_upper_triangular(&v);
            return Some((v, vinv));
        }

        // Repeated eigenvalues: check diagonalizability via rank condition.
        // For each distinct eigenvalue λ with multiplicity m, geometric
        // multiplicity is n - rank(T - λI): diagonalizable iff it equals m,
        // i.e. rank(T - λI) == n - m. Anything else (including
        // rank > n - m, the defective case) rejects the decomposition.
        let n = t.dim;
        let eigvals = Self::eigenvalues_from_triangular(t);

        // Group eigenvalues by proximity
        let mut processed = vec![false; n];
        for (i, eig_i) in eigvals.iter().enumerate() {
            if processed[i] {
                continue;
            }
            // Find all eigenvalues close to eig_i
            let mut multiplicity = 0_usize;
            for j in i..n {
                if eig_i.sub(&eigvals[j]).abs().total_cmp(&tol) == Ordering::Less {
                    processed[j] = true;
                    multiplicity += 1;
                }
            }

            // Compute rank(T - λ_i * I)
            let mut shifted = t.clone();
            for k in 0..n {
                let val = shifted.get(k, k).sub(eig_i);
                shifted.set(k, k, val);
            }
            let scalar_ratio = rank_upper_triangular(&shifted, &tol);
            if scalar_ratio != n - multiplicity {
                return None; // defective
            }
        }

        // Matrix is diagonalizable with repeated eigenvalues.
        // For each distinct eigenvalue, find nullspace basis of (T - λI).
        let mut v = MatC::identity(n);
        processed = vec![false; n];
        for (i, eig_i) in eigvals.iter().enumerate() {
            if processed[i] {
                continue;
            }
            // Find cluster: all j where λ_j ≈ λ_i
            let mut cluster: Vec<usize> = Vec::new();
            for j in 0..n {
                if eig_i.sub(&eigvals[j]).abs().total_cmp(&tol) == Ordering::Less && !processed[j] {
                    cluster.push(j);
                    processed[j] = true;
                }
            }
            let m = cluster.len();
            if m == 1 {
                // Simple eigenvalue — already handled by identity initialization
                continue;
            }

            // Solve for nullspace basis vectors of (T - λI).
            // The nullspace has dimension m. For each cluster position idx,
            // set v[idx][idx] = 1 and solve for earlier components in the cluster.
            let lambda = eig_i;
            for &idx in &cluster {
                for &j in &cluster {
                    if j < idx {
                        let mut sum = Cmplx::zero();
                        for &l in &cluster {
                            if l > j && l <= idx {
                                sum += &t.get(j, l).mul(v.get(l, idx));
                            }
                        }
                        let denom = t.get(j, j).sub(lambda);
                        if denom.abs().total_cmp(&tol) == Ordering::Greater {
                            v.set(j, idx, sum.neg().div(&denom));
                        }
                    }
                }
            }
        }

        let vinv = inverse_upper_triangular(&v);
        Some((v, vinv))
    }

    /// Solve `M·X = I` via LU with partial pivoting, or `None` if singular.
    ///
    /// The standard dense inverse: backward-stable, eigenvector-free, and
    /// exact for small integer matrices. Singularity is scale-aware
    /// (pivot `<= eps * max|M|` declines). This is the inversion path;
    /// spectral evaluation stays on Schur/Parlett.
    pub fn lu_solve(mat: &MatC) -> Option<MatC> {
        let n = mat.dim;
        if n == 0 {
            return None;
        }
        let mut lu = mat.clone();
        let mut perm: Vec<usize> = (0..n).collect();
        let mut max_abs = Real::zero();
        for i in 0..n {
            for j in 0..n {
                let a = lu.get(i, j).abs();
                if a.total_cmp(&max_abs) == Ordering::Greater {
                    max_abs = a;
                }
            }
        }
        let tol = Real::epsilon() * &max_abs;
        for k in 0..n {
            let mut piv = k;
            let mut piv_abs = lu.get(k, k).abs();
            for i in k + 1..n {
                let a = lu.get(i, k).abs();
                if a.total_cmp(&piv_abs) == Ordering::Greater {
                    piv = i;
                    piv_abs = a;
                }
            }
            if piv_abs.total_cmp(&tol) != Ordering::Greater {
                return None;
            }
            if piv != k {
                for j in 0..n {
                    let tmp = lu.get(k, j).clone();
                    lu.set(k, j, lu.get(piv, j).clone());
                    lu.set(piv, j, tmp);
                }
                perm.swap(k, piv);
            }
            for i in k + 1..n {
                let factor = lu.get(i, k).div(lu.get(k, k));
                lu.set(i, k, factor.clone());
                for j in k + 1..n {
                    let sub = factor.mul(lu.get(k, j));
                    lu.set(i, j, lu.get(i, j).sub(&sub));
                }
            }
        }
        // Back/forward substitution per column: rhs = P·e_col.
        let mut inv = MatC::zeros(n);
        for col in 0..n {
            let mut y: Vec<Cmplx> = Vec::with_capacity(n);
            for (i, p) in perm.iter().enumerate() {
                let mut acc = if *p == col {
                    Cmplx::one()
                } else {
                    Cmplx::zero()
                };
                for (j, yj) in y.iter().enumerate() {
                    acc = acc.sub(&lu.get(i, j).mul(yj));
                }
                y.push(acc);
            }
            // Back-substitution writes rows n-1 down to 0, so already-solved
            // rows are readable directly from `inv`.
            for i in (0..n).rev() {
                let mut acc = y[i].clone();
                for j in i + 1..n {
                    acc = acc.sub(&lu.get(i, j).mul(inv.get(j, col)));
                }
                inv.set(i, col, acc.div(lu.get(i, i)));
            }
        }
        Some(inv)
    }
}

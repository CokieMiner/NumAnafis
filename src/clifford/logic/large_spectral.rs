//! Spectral function application for large Clifford algebras.
//!
//! Block-Parlett evaluation on the Schur form with eigenvalue clustering,
//! plus the public `apply_*` entry points used by the small-algebra
//! dispatcher on fallback.

use core::cmp::Ordering;

use alloc::vec::Vec;

use crate::{
    traits::Numeric,
    types::{IntType, Real},
};

use super::{CliffordNumber, Cmplx, EigenSolver, Jet, LargeEmbed, MatC, SpectralFn};

// Heuristic native work limits; acceptance uses the terminal-window test below.
// TODO(mp): derive work limits from backend precision, conditioning, and the
// requested error; validate terminal-window acceptance across precisions.
const TAYLOR_EXTRA_TERMS: usize = 64;
const TAYLOR_TAIL_EXTRA_TERMS: usize = 8;

/// Heuristic eigenvalue-clustering threshold: sqrt(eps) times the largest entry magnitude.
/// Schur reordering groups connected components of the resulting proximity graph.
fn cluster_tol(t: &MatC) -> Real {
    let mut scale = Real::zero();
    for i in 0..t.dim {
        for j in 0..t.dim {
            let a = t.get(i, j).abs();
            if a.total_cmp(&scale) == Ordering::Greater {
                scale = a;
            }
        }
    }
    Real::epsilon().sqrt() * scale
}

/// Write a computed diagonal block into the result matrix.
fn write_block(f: &mut MatC, f_block: &MatC, start: usize) {
    let n = f_block.dim;
    for i in 0..n {
        for j in i..n {
            f.set(start + i, start + j, f_block.get(i, j).clone());
        }
    }
}

/// Exact-zero test for a matrix (nilpotent termination, no tolerance).
fn is_zero_mat(m: &MatC) -> bool {
    for i in 0..m.dim {
        for j in 0..m.dim {
            if !m.get(i, j).is_zero() {
                return false;
            }
        }
    }
    true
}

/// Maximum entry magnitude of a matrix.
fn max_abs_mat(m: &MatC) -> Real {
    let mut best = Real::zero();
    for i in 0..m.dim {
        for j in 0..m.dim {
            let a = m.get(i, j).abs();
            if a.total_cmp(&best) == Ordering::Greater {
                best = a;
            }
        }
    }
    best
}

/// Evaluate an AST on an upper triangular Schur matrix via Block Parlett recurrence.
pub fn apply_ast_schur(t: &MatC, ast: &SpectralFn, blocks: &[(usize, usize)]) -> Option<MatC> {
    let n = t.dim;
    let mut f = MatC::zeros(n);
    let tol = cluster_tol(t);

    // Evaluate the Taylor series for each diagonal block.
    for &(start, end) in blocks {
        let block_size = end - start + 1;
        let lambda = t.get(start, start);

        // The Taylor remainder N = T_block - lambda I includes diagonal spread.
        let mut n_mat = MatC::zeros(block_size);
        let mut n_zero = true;
        for i in 0..block_size {
            for j in i..block_size {
                let mut v = t.get(start + i, start + j).clone();
                if i == j {
                    v = v.sub(lambda);
                }
                if !v.is_zero() {
                    n_zero = false;
                }
                n_mat.set(i, j, v);
            }
        }

        let constant = ast.eval(lambda);
        if !constant.0.is_finite() || !constant.1.is_finite() {
            return None;
        }
        let mut f_block = MatC::zeros(block_size);
        for i in 0..block_size {
            f_block.set(i, i, constant.clone());
        }
        if n_zero {
            write_block(&mut f, &f_block, start);
            continue;
        }
        // Strictly upper triangular remainders terminate exactly at block_size.
        let nilpotent = (0..block_size).all(|i| n_mat.get(i, i).is_zero());
        let degree = if nilpotent {
            block_size
        } else {
            block_size.checked_add(TAYLOR_EXTRA_TERMS)?
        };
        let result_jet = ast.eval_jet(&Jet::variable(lambda, degree));
        // General case: accumulate f(λ)I + c_1 N + c_2 N² + ….
        // Terms below `block_size` are always taken (nilpotent-exact region:
        // skipping a zero early term would miss later ones, e.g. f(z) = z²).
        // Evaluate the full work budget and test the terminal window below.
        let taylor_tol = Real::epsilon() * Real::from_int(100);
        let mut n_power = MatC::identity(block_size);
        let mut terminated = false;
        let mut tail = Real::zero();
        let window = block_size.checked_add(TAYLOR_TAIL_EXTRA_TERMS)?;
        for (m, c) in result_jet.coeffs.iter().enumerate().skip(1) {
            n_power = n_power.mul(&n_mat);
            if is_zero_mat(&n_power) {
                terminated = true;
                break;
            }
            if !c.0.is_finite() || !c.1.is_finite() {
                return None;
            }
            let term = n_power.scale(c);
            if !term.data.iter().all(|v| v.0.is_finite() && v.1.is_finite()) {
                return None;
            }
            f_block = f_block.add(&term);
            if m + window >= result_jet.coeffs.len() {
                tail = tail + max_abs_mat(&term);
            }
        }
        // A full terminal window guards against isolated zero Taylor terms.
        // This is a numerical acceptance criterion, not a certified error bound.
        let scale = max_abs_mat(&f_block) + constant.abs();
        if !terminated && (nilpotent || tail > scale * taylor_tol) {
            return None;
        }

        write_block(&mut f, &f_block, start);
    }

    // Apply the Parlett recurrence to off-diagonal blocks.
    for j in 1..n {
        for i in (0..j).rev() {
            // Check if i and j are in the same block
            let same_block = blocks.iter().any(|&(scalar, e)| i >= scalar && j <= e);
            if same_block {
                continue;
            }

            let mut sum = Cmplx::zero();
            for k in i + 1..j {
                let term1 = f.get(i, k).mul(t.get(k, j));
                let term2 = t.get(i, k).mul(f.get(k, j));
                sum += &term1.sub(&term2);
            }
            sum += &f.get(i, i).mul(t.get(i, j));
            sum = sum.sub(&t.get(i, j).mul(f.get(j, j)));

            let denom = t.get(i, i).sub(t.get(j, j));
            if denom.abs().total_cmp(&tol) == Ordering::Less {
                // Coincident eigenvalues across distinct blocks invalidate this recurrence.
                return None;
            }

            f.set(i, j, sum.div(&denom));
        }
    }

    Some(f)
}

/// Reorders the Schur decomposition so that clustered eigenvalues are
/// contiguous on the diagonal (Bai–Demmel swaps via Givens rotations).
///
/// Clustering is transitive (union-find) under the shared [`cluster_tol`],
/// so chains of nearby eigenvalues stay together exactly as the Parlett
/// stage groups them.
pub fn reorder_schur(t: &mut MatC, q: &mut MatC) -> Vec<(usize, usize)> {
    #[expect(
        clippy::missing_const_for_fn,
        reason = "Slice indexing and mutation are not const-compatible"
    )]
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let n = t.dim;
    let tol = cluster_tol(t);
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in i + 1..n {
            let close = t.get(i, i).sub(t.get(j, j)).abs().total_cmp(&tol) != Ordering::Greater;
            if close {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }
    // Canonical cluster IDs in order of first appearance.
    let mut cluster_id = alloc::vec![0_usize; n];
    let mut root_to_id: Vec<(usize, usize)> = Vec::new();
    for (i, slot) in cluster_id.iter_mut().enumerate() {
        let scalar_ratio = find(&mut parent, i);
        let id = if let Some((_, id)) = root_to_id.iter().find(|(root, _)| *root == scalar_ratio) {
            *id
        } else {
            let id = root_to_id.len() + 1;
            root_to_id.push((scalar_ratio, id));
            id
        };
        *slot = id;
    }

    // Sort the eigenvalues so that elements of the same cluster are contiguous.
    let mut swapped = true;
    while swapped {
        swapped = false;
        for k in 0..n.saturating_sub(1) {
            if cluster_id[k] > cluster_id[k + 1] {
                swap_adjacent_schur(t, q, k);
                cluster_id.swap(k, k + 1);
                swapped = true;
            }
        }
    }
    let mut blocks = Vec::new();
    let mut start = 0;
    for i in 1..n {
        if cluster_id[i] != cluster_id[i - 1] {
            blocks.push((start, i - 1));
            start = i;
        }
    }
    if n != 0 {
        blocks.push((start, n - 1));
    }
    blocks
}

#[expect(
    clippy::many_single_char_names,
    reason = "Matrix math conventionally uses single letters"
)]
fn swap_adjacent_schur(t: &mut MatC, q: &mut MatC, k: usize) {
    let t11 = t.get(k, k).clone();
    let t12 = t.get(k, k + 1).clone();
    let t22 = t.get(k + 1, k + 1).clone();

    let x = t12;
    let y = t22.sub(&t11);

    let norm_sqr = x.abs_sq() + y.abs_sq();
    if norm_sqr.total_cmp(&Real::zero()) == Ordering::Equal {
        return; // Already decoupled
    }

    let norm = norm_sqr.sqrt();
    let norm_cmplx = Cmplx(norm, Real::zero());
    let c = x.div(&norm_cmplx);
    let s_rot = y.div(&norm_cmplx);
    let c_conj = c.conj();
    let s_conj = s_rot.conj();

    let n = t.dim;

    // T = T * G
    for i in 0..n {
        let col_k = t.get(i, k).clone();
        let col_k1 = t.get(i, k + 1).clone();

        let new_col_k = col_k.mul(&c).add(&col_k1.mul(&s_rot));
        let new_col_k1 = col_k.mul(&s_conj.neg()).add(&col_k1.mul(&c_conj));

        t.set(i, k, new_col_k);
        t.set(i, k + 1, new_col_k1);
    }

    // T = G^* * T
    for j in 0..n {
        let row_k = t.get(k, j).clone();
        let row_k1 = t.get(k + 1, j).clone();

        let new_row_k = c_conj.mul(&row_k).add(&s_conj.mul(&row_k1));
        let new_row_k1 = s_rot.neg().mul(&row_k).add(&c.mul(&row_k1));

        t.set(k, j, new_row_k);
        t.set(k + 1, j, new_row_k1);
    }
    t.set(k + 1, k, Cmplx::zero());

    // Q = Q * G
    for i in 0..n {
        let col_k = q.get(i, k).clone();
        let col_k1 = q.get(i, k + 1).clone();

        let new_col_k = col_k.mul(&c).add(&col_k1.mul(&s_rot));
        let new_col_k1 = col_k.mul(&s_conj.neg()).add(&col_k1.mul(&c_conj));

        q.set(i, k, new_col_k);
        q.set(i, k + 1, new_col_k1);
    }
}

/// Check `T·V ≈ V·D` for a triangular `T`, its purported eigenvectors `V`,
/// and `D = diag(T)`, against a scale-aware tolerance.
///
/// Ill-conditioned eigenvector matrices pass the rank test yet amplify
/// rounding in the `f(D)` sandwich by `‖V‖`; the residual below scales
/// with that amplification, so garbage factors are declined here rather
/// than silently consumed.
fn eigen_decomp_valid(t: &MatC, v: &MatC) -> bool {
    let n = t.dim;
    let mut scale = Real::zero();
    for i in 0..n {
        for j in 0..n {
            let a = t.get(i, j).abs();
            if a.total_cmp(&scale) == Ordering::Greater {
                scale = a;
            }
        }
    }
    let nf =
        Real::from_int(IntType::try_from(n).expect("matrix dimension fits the integer backend"));
    let tol = Real::epsilon() * &scale * &nf * &nf * &Real::from_int(100);
    for i in 0..n {
        for j in 0..n {
            // (T·V)[i][j]
            let mut tv = Cmplx::zero();
            for k in 0..n {
                tv = tv.add(&t.get(i, k).mul(v.get(k, j)));
            }
            // (V·D)[i][j] = V[i][j] * T[j][j]
            let vd = v.get(i, j).mul(t.get(j, j));
            if tv.sub(&vd).abs().total_cmp(&tol) == Ordering::Greater {
                return false;
            }
        }
    }
    true
}

/// Verify a projected multivector re-embeds to its source matrix.
///
/// Spectral images need not lie in the represented algebra; the trace
/// projection then returns coefficients that do not satisfy the function
/// equation. Callers decline (NaN) when the residual exceeds tolerance.
fn projection_valid(mat: &MatC, mv: &CliffordNumber, basis: &[MatC]) -> bool {
    let rt = LargeEmbed::to_matrix(mv, basis);
    if rt.dim != mat.dim {
        return false;
    }
    let mut scale = Real::zero();
    for i in 0..mat.dim {
        for j in 0..mat.dim {
            let a = mat.get(i, j).abs();
            if a.total_cmp(&scale) == Ordering::Greater {
                scale = a;
            }
        }
    }
    let tol = Real::epsilon() * scale * Real::from_int(100);
    for i in 0..mat.dim {
        for j in 0..mat.dim {
            if mat.get(i, j).sub(rt.get(i, j)).abs().total_cmp(&tol) == Ordering::Greater {
                return false;
            }
        }
    }
    true
}
// LargeSpectral — public large-algebra application entries

/// Large-algebra spectral application: the public entry points that
/// evaluate analytic functions and ASTs on the Schur form.
pub struct LargeSpectral;

impl LargeSpectral {
    /// Apply an analytic function `f: Cmplx → Cmplx` via Schur decomposition.
    ///
    /// Returns NaN if embedding, decomposition, or representability fails.
    pub fn large_apply_via_spectral_closure<F>(mv: &CliffordNumber, f: &F) -> CliffordNumber
    where
        F: Fn(&Cmplx) -> Cmplx,
    {
        let gens = mv.generator_set();
        let Some(basis) = LargeEmbed::blade_basis(gens) else {
            return CliffordNumber::nan(gens.clone());
        };
        let mat = LargeEmbed::to_matrix(mv, &basis);

        let Some((q, t)) = EigenSolver::schur_decomposition(&mat) else {
            return CliffordNumber::nan(gens.clone());
        };

        let dim = t.dim;

        // Attempt full eigendecomposition from the Schur form.
        // The factors must satisfy T·V ≈ V·D, or the f(D) sandwich below
        // silently produces garbage: decline instead.
        let Some((v, vinv)) = EigenSolver::eigendecompose_from_schur(&t) else {
            return CliffordNumber::nan(gens.clone());
        };
        if !eigen_decomp_valid(&t, &v) {
            return CliffordNumber::nan(gens.clone());
        }

        // Build f(D) — diagonal matrix of f(eigenvalues)
        let eigvals = EigenSolver::eigenvalues_from_triangular(&t);
        let mut fd_mat = MatC::zeros(dim);
        for (i, eig) in eigvals.iter().enumerate() {
            let f_eig = f(eig);
            if !f_eig.0.is_finite() || !f_eig.1.is_finite() {
                return CliffordNumber::nan(gens.clone());
            }
            fd_mat.set(i, i, f_eig);
        }

        // result = Q * V * f(D) * V^{-1} * Q^H
        let qv = q.mul(&v);
        let fd_vinv = fd_mat.mul(&vinv);
        let qv_fd_vinv = qv.mul(&fd_vinv);
        let result_mat = qv_fd_vinv.mul(&q.dagger());

        let projected = LargeEmbed::from_matrix(&result_mat, gens, &basis);
        if projection_valid(&result_mat, &projected, &basis) {
            return projected;
        }
        CliffordNumber::nan(gens.clone())
    }

    pub fn large_apply_via_spectral(mv: &CliffordNumber, ast: &SpectralFn) -> CliffordNumber {
        let gens = mv.generator_set();
        let Some(basis) = LargeEmbed::blade_basis(gens) else {
            return CliffordNumber::nan(gens.clone());
        };
        let mat = LargeEmbed::to_matrix(mv, &basis);

        let Some((mut q, mut t)) = EigenSolver::schur_decomposition(&mat) else {
            return CliffordNumber::nan(gens.clone());
        };

        let blocks = reorder_schur(&mut t, &mut q);

        let Some(f_mat) = apply_ast_schur(&t, ast, &blocks) else {
            return CliffordNumber::nan(gens.clone());
        };

        // result = Q * F * Q^H
        let q_f = q.mul(&f_mat);
        let result_mat = q_f.mul(&q.dagger());

        let projected = LargeEmbed::from_matrix(&result_mat, gens, &basis);
        if projection_valid(&result_mat, &projected, &basis) {
            return projected;
        }
        CliffordNumber::nan(gens.clone())
    }
}

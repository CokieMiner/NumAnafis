//! Complex-matrix Clifford embeddings with signature-dependent dimensions.
//!
//! Dense complex matrix type with Kronecker-built gamma matrices and
//! multivector projection. Eigensolving lives in [`super::large_eigen`],
//! spectral application in [`super::large_spectral`].

use alloc::{vec, vec::Vec};

use crate::types::{RationalType, Real};

use super::{BladeAlgebra, CliffordNumber, Cmplx, GeneratorSet};
// Custom Complex Matrix

#[derive(Clone, Debug)]
pub struct MatC {
    pub dim: usize,
    pub data: Vec<Cmplx>,
}

impl MatC {
    pub fn zeros(dim: usize) -> Self {
        let mut data = Vec::with_capacity(dim * dim);
        let z = Cmplx::zero();
        for _ in 0..dim * dim {
            data.push(z.clone());
        }
        Self { dim, data }
    }

    pub fn identity(dim: usize) -> Self {
        let mut mat = Self::zeros(dim);
        let one = Cmplx::one();
        for i in 0..dim {
            mat.set(i, i, one.clone());
        }
        mat
    }

    pub fn get(&self, row: usize, col: usize) -> &Cmplx {
        &self.data[row * self.dim + col]
    }

    pub fn set(&mut self, row: usize, col: usize, val: Cmplx) {
        self.data[row * self.dim + col] = val;
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = Self::zeros(self.dim);
        for i in 0..self.data.len() {
            result.data[i] = self.data[i].add(&other.data[i]);
        }
        result
    }

    pub fn mul(&self, other: &Self) -> Self {
        let dim = self.dim;
        let mut result = Self::zeros(dim);
        for i in 0..dim {
            for j in 0..dim {
                let mut sum = Cmplx::zero();
                for k in 0..dim {
                    sum += &self.get(i, k).mul(other.get(k, j));
                }
                result.set(i, j, sum);
            }
        }
        result
    }

    pub fn scale(&self, scalar: &Cmplx) -> Self {
        let mut result = Self::zeros(self.dim);
        for i in 0..self.data.len() {
            result.data[i] = self.data[i].mul(scalar);
        }
        result
    }

    pub fn dagger(&self) -> Self {
        let mut result = Self::zeros(self.dim);
        for i in 0..self.dim {
            for j in 0..self.dim {
                result.set(j, i, self.get(i, j).conj());
            }
        }
        result
    }
}
// Kronecker product

fn kron(a: &MatC, b: &MatC) -> MatC {
    let dim = a.dim * b.dim;
    let mut result = MatC::zeros(dim);
    for i in 0..a.dim {
        for j in 0..a.dim {
            let a_ij = a.get(i, j);
            if a_ij.is_zero() {
                continue;
            }
            let base_row = i * b.dim;
            let base_col = j * b.dim;
            for p in 0..b.dim {
                for q in 0..b.dim {
                    let b_pq = b.get(p, q);
                    if b_pq.is_zero() {
                        continue;
                    }
                    result.set(base_row + p, base_col + q, a_ij.mul(b_pq));
                }
            }
        }
    }
    result
}

fn kron_list(factors: &[MatC]) -> MatC {
    if factors.is_empty() {
        let mut m = MatC::zeros(1);
        m.set(0, 0, Cmplx::one());
        return m;
    }
    let mut result = factors[0].clone();
    for f in &factors[1..] {
        result = kron(&result, f);
    }
    result
}
// Pauli matrices

fn pauli_sx() -> MatC {
    let mut m = MatC::zeros(2);
    m.set(0, 1, Cmplx::one());
    m.set(1, 0, Cmplx::one());
    m
}

fn pauli_sy() -> MatC {
    let mut m = MatC::zeros(2);
    let i = Cmplx::i();
    let neg_i = i.neg();
    m.set(0, 1, neg_i);
    m.set(1, 0, i);
    m
}

fn pauli_sz() -> MatC {
    let mut m = MatC::zeros(2);
    m.set(0, 0, Cmplx::one());
    m.set(1, 1, Cmplx::one().neg());
    m
}

fn identity2() -> MatC {
    MatC::identity(2)
}
// Euclidean gamma matrices (Weyl / chiral representation)

fn euclidean_gammas(n: usize) -> Vec<MatC> {
    let k = n.div_ceil(2);
    let sx = pauli_sx();
    let sy = pauli_sy();
    let sz = pauli_sz();
    let i2 = identity2();

    let mut gammas = Vec::with_capacity(n);
    for pair in 0..n.div_euclid(2) {
        let mut x_factors: Vec<MatC> = Vec::with_capacity(k);
        let mut y_factors: Vec<MatC> = Vec::with_capacity(k);
        for _ in 0..pair {
            x_factors.push(sz.clone());
            y_factors.push(sz.clone());
        }
        x_factors.push(sx.clone());
        y_factors.push(sy.clone());
        for _ in 0..(k - pair - 1) {
            x_factors.push(i2.clone());
            y_factors.push(i2.clone());
        }
        gammas.push(kron_list(&x_factors));
        gammas.push(kron_list(&y_factors));
    }
    if n % 2 == 1 {
        let factors = vec![sz; k];
        gammas.push(kron_list(&factors));
    }
    gammas
}
// Metric-adjusted generator matrices

fn generator_matrices(gens: &GeneratorSet) -> Option<Vec<MatC>> {
    // J=P*g_last=i*I. Since P^2=s, g_last=i*P/s.
    // This selects the same complex orientation as public scalar promotion.
    if gens.complex_blade().is_some() {
        let n = gens.len();
        let prefix = GeneratorSet::from_validated(
            (0..n - 1)
                .map(|index| (gens.id_at(index), gens.metric_at(index)))
                .collect(),
        );
        let mut result = if n == 1 {
            Vec::new()
        } else {
            generator_matrices(&prefix)?
        };
        let mut product = MatC::identity(result.first().map_or(1, |matrix| matrix.dim));
        for matrix in &result {
            product = product.mul(matrix);
        }
        let blade = (1_usize << (n - 1)) - 1;
        let (_, sign) = BladeAlgebra::mul_blades(&prefix, blade, blade)?;
        let factor = if sign == 1 {
            Cmplx::i()
        } else {
            Cmplx::i().neg()
        };
        result.push(product.scale(&factor));
        return Some(result);
    }

    let mut num_nilpotent = 0;
    for i in 0..gens.len() {
        if gens.metric_at(i) == 0 {
            num_nilpotent += 1;
        }
    }

    let euclid = euclidean_gammas(gens.len() + num_nilpotent);
    let mut result = Vec::with_capacity(gens.len());

    let mut gamma_idx = 0;
    for i in 0..gens.len() {
        match gens.metric_at(i) {
            1 => {
                result.push(euclid[gamma_idx].clone());
                gamma_idx += 1;
            }
            -1 => {
                result.push(euclid[gamma_idx].scale(&Cmplx::i()));
                gamma_idx += 1;
            }
            0 => {
                // n = (gamma_1 + i * gamma_2) / 2
                // Since gamma_1 and gamma_2 anticommute and square to 1, n^2 = 0
                let g1 = &euclid[gamma_idx];
                let g2 = &euclid[gamma_idx + 1];
                let i_g2 = g2.scale(&Cmplx::i());
                let mut n = g1.add(&i_g2);
                let half = Cmplx::new(
                    Real::from_rational(RationalType::from_parts(1, 2)),
                    Real::zero(),
                );
                n = n.scale(&half);
                result.push(n);
                gamma_idx += 2;
            }
            _ => return None,
        }
    }
    Some(result)
}
// LargeEmbed — MatC embedding for n >= 4 generators

/// Large-algebra matrix embedding: gamma-matrix blade basis and
/// multivector projection to and from dense complex matrices.
pub struct LargeEmbed;

impl LargeEmbed {
    pub fn blade_basis(gens: &GeneratorSet) -> Option<Vec<MatC>> {
        let gen_mats = generator_matrices(gens)?;
        // Zero-generator algebras use scalar dispatch instead of this embedding.
        let first = gen_mats.first()?;
        let n = gen_mats.len();
        let count = 1_usize << n;
        let dim = first.dim;
        let mut basis = Vec::with_capacity(count);
        basis.push(MatC::identity(dim));
        for i in 1..count {
            let mut blade = MatC::identity(dim);
            for (g, gm) in gen_mats.iter().enumerate() {
                if (i >> g) & 1 == 1 {
                    blade = blade.mul(gm);
                }
            }
            basis.push(blade);
        }
        Some(basis)
    }
    // CliffordNumber ↔ MatC

    pub fn to_matrix(mv: &CliffordNumber, basis: &[MatC]) -> MatC {
        let dim = basis[0].dim;
        let mut mat = MatC::zeros(dim);
        for (blade, bm) in basis.iter().enumerate() {
            let coeff = mv.coeff(blade);
            if coeff.is_zero() {
                continue;
            }
            for (target, source) in mat.data.iter_mut().zip(&bm.data) {
                *target += &source.scale(coeff);
            }
        }
        mat
    }

    pub fn from_matrix(mat: &MatC, gens: &GeneratorSet, basis: &[MatC]) -> CliffordNumber {
        let n = gens.len();
        let count = match u8::try_from(n) {
            Ok(n_u8) => BladeAlgebra::coeff_count(n_u8).unwrap_or(0),
            Err(_) => return CliffordNumber::nan(gens.clone()),
        };
        let mut mv = CliffordNumber::zero_unchecked(gens.clone());
        for (blade, bm) in basis.iter().enumerate().take(count) {
            // Re tr(A B*) = sum(Re(Aij)*Re(Bij)+Im(Aij)*Im(Bij)).
            // Only the trace is needed: no dense matrix products or temporaries.
            let mut numer = Real::zero();
            let mut denom = Real::zero();
            for (a, b) in mat.data.iter().zip(&bm.data) {
                numer = numer + &a.0 * &b.0 + &a.1 * &b.1;
                denom = denom + b.abs_sq();
            }

            let c = if denom.is_zero() {
                Real::zero()
            } else {
                numer / denom
            };
            if !c.is_zero() {
                mv.coeffs_mut_slice()[blade] = c;
            }
        }
        mv
    }
}

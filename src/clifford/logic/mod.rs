//! Multivector logic: dense types, CGA constructors, arithmetic, matrix
//! embeddings, and spectral (functional-calculus) evaluation.

mod arithmetic;
mod clifford_numeric;
mod constructors;
mod core;
mod core_products;
mod display;
mod fast;
mod gen_id;
mod large_eigen;
mod large_mat;
mod large_spectral;
mod matrix;
mod spectral;
mod spectral_ast;
mod spectral_jet;
mod types;

pub use constructors::Cga;
pub use fast::FastClifford;
pub use gen_id::{E_MINUS, E_PLUS, E1, E2, E3};
pub use large_eigen::EigenSolver;
pub use large_mat::{LargeEmbed, MatC};
pub use large_spectral::LargeSpectral;
pub use matrix::{Cmplx, Mat2C, SmallEmbed};
pub use spectral::SpectralDispatch;
pub use spectral_ast::SpectralFn;
pub use spectral_jet::Jet;
pub use types::{
    BladeAlgebra, CliffordCoeffs, CliffordNumber, GeneratorSet, INLINE_COEFF_COUNT,
    INLINE_GENERATOR_LIMIT,
};

#[cfg(test)]
mod tests;

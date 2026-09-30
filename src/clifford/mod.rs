//! Dense multivectors for Clifford (geometric) algebras with labeled generators.
//!
//! This module implements **Conformal Geometric Algebra** CGA Cl(4,1) with 5
//! base generators: `e1, e2, e3` (Euclidean 3-space), `e+` (extra positive),
//! and `e−` (extra negative).  All derived objects (quaternions, complex unit,
//! conformal origin/infinity) are constructed as products of these generators.

#![expect(
    clippy::indexing_slicing,
    reason = "Blade indices and matrix dimensions are derived from validated signatures and coefficient lengths"
)]

mod api;
mod logic;

pub use api::{Cga, CliffordNumber, FastClifford, GeneratorSet};

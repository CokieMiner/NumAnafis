//! Stable IDs for the CGA Cl(4,1) base generators and user-defined extensions.
//!
//! The five fundamental generators are ordered by ID so that a
//! [`crate::clifford::GeneratorSet`] built from all five is already sorted.
//!
//! IDs 0–255 are reserved by the library.
//! User-defined generators should start at 256.

/// Euclidean basis vector `e1` — `e1² = +1`.
pub const E1: u32 = 0;
/// Euclidean basis vector `e2` — `e2² = +1`.
pub const E2: u32 = 1;
/// Euclidean basis vector `e3` — `e3² = +1`.
pub const E3: u32 = 2;
/// Extra positive dimension `e+` — `e+² = +1`.
pub const E_PLUS: u32 = 3;
/// Extra negative dimension `e−` — `e−² = −1`.
pub const E_MINUS: u32 = 4;

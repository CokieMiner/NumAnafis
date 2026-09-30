//! Numeric values for symbolic evaluation, with exact native integer/rational
//! arithmetic, real-component complex values, and optional Clifford algebras.
//!
//! [`Number`] stores every supported domain in one evaluator value. [`Real`]
//! supplies its integer, rational, and approximate tiers. [`Numeric`] registers
//! supported elementary functions; special kernels remain available separately
//! through [`F32Ext`] and [`F64Ext`].

#![cfg_attr(not(feature = "std"), no_std)]

extern crate alloc;
extern crate libm;

mod error;
mod ext;
mod special;
mod traits;
mod types;

#[cfg(feature = "clifford")]
mod clifford;

#[cfg(feature = "python")]
mod python_binding;

// Primary re-exports at root
pub use error::NumAnafisError;
pub use ext::{F32Ext, F64Ext};
pub use traits::Numeric;
pub use types::{Complex, FloatType, IntType, Number, RationalType, Real, c, i, n, r};

#[cfg(feature = "clifford")]
pub use clifford::{Cga, CliffordNumber, FastClifford, GeneratorSet};

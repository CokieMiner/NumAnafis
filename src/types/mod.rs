//! Native real, complex, and optional Clifford numeric representations.

mod complex;
mod complex_numeric;
mod constructors;
mod float;
mod int;
mod number;
mod number_inspect;
mod number_numeric;
mod number_ops;
mod rational;
mod real;
mod real_compare;
mod real_numeric;
mod real_ops;

pub use complex::Complex;
pub use constructors::{c, i, n, r};
pub use float::{FloatMath, FloatType};
pub use int::{IntMath, IntType};
pub use number::{Number, NumberRepr};
pub use rational::{RationalMath, RationalType};
pub use real::Real;

#[cfg(test)]
mod tests;

//! Error types for `num-anafis`.

#[cfg(feature = "std")]
use core::error::Error;
use core::fmt::{Display, Formatter, Result};

use alloc::boxed::Box;

/// Error type for invalid `num-anafis` operations and constructors.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum NumAnafisError {
    /// A real-only representation was required.
    ExpectedReal,
    /// The algebra has no declared central square-minus-one pseudoscalar.
    UnsupportedComplexStructure,
    /// An ordinary real or complex scalar was required.
    ExpectedScalar,
    /// Active generators exceed the signature capacity.
    ActiveGeneratorsExceedSignature(Box<SignatureMismatchError>),
    /// Generator index is out of bounds for the active algebra.
    GeneratorIndexOutOfRange(Box<IndexOutOfRangeError>),
    /// A generator label appears more than once.
    DuplicateGenerator {
        /// Repeated generator label.
        id: u32,
    },
    /// A diagonal metric is outside `{-1, 0, 1}`.
    InvalidGeneratorMetric {
        /// Generator label.
        id: u32,
        /// Invalid metric.
        metric: i8,
    },
    /// Dense coefficients length does not match `2^n`.
    DenseCoefficientLengthMismatch(Box<DenseLengthError>),
    /// A sparse term names a blade outside the active algebra.
    BladeIndexOutOfRange {
        /// Invalid blade mask.
        blade: usize,
        /// Active blade count.
        count: usize,
    },
    /// Active generators exceed platform indexable width for `usize` bitmasks.
    ActiveGeneratorsTooLargeForPlatform {
        /// Requested active generator count.
        active: u8,
    },
    /// Attempted to cast a `CliffordNumber` from an incompatible algebra.
    MismatchedGeneratorSet,
}

/// Error indicating active generators exceed signature available generators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct SignatureMismatchError {
    /// Active generator count.
    pub active: u8,
    /// Available generator count.
    pub available: u8,
}

/// Error indicating generator index is out of valid range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct IndexOutOfRangeError {
    /// The out-of-range index.
    pub index: u8,
    /// Active generator count.
    pub active: u8,
}

/// Error details for length mismatches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct DenseLengthError {
    /// Expected coefficient count.
    pub expected: usize,
    /// Provided coefficient count.
    pub found: usize,
}

impl Display for NumAnafisError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        match *self {
            Self::UnsupportedComplexStructure => {
                write!(f, "the algebra has no compatible central complex structure")
            }
            Self::ExpectedReal => write!(f, "a real component is required"),
            Self::ExpectedScalar => write!(f, "a scalar coefficient is required"),
            Self::ActiveGeneratorsExceedSignature(ref e) => {
                let (active, available) = (e.active, e.available);
                write!(
                    f,
                    "active generators ({active}) exceed signature generators ({available})"
                )
            }
            Self::GeneratorIndexOutOfRange(ref e) => {
                let (index, active) = (e.index, e.active);
                write!(f, "generator index {index} out of range for n={active}")
            }
            Self::DuplicateGenerator { id } => {
                write!(f, "duplicate generator label {id}")
            }
            Self::InvalidGeneratorMetric { id, metric } => {
                write!(f, "generator {id} has invalid diagonal metric {metric}")
            }
            Self::DenseCoefficientLengthMismatch(ref e) => {
                let (expected, found) = (e.expected, e.found);
                write!(
                    f,
                    "dense coefficient length mismatch: expected {expected}, got {found}"
                )
            }
            Self::BladeIndexOutOfRange { blade, count } => {
                write!(
                    f,
                    "blade mask {blade} is outside the active range 0..{count}"
                )
            }
            Self::ActiveGeneratorsTooLargeForPlatform { active } => {
                write!(
                    f,
                    "active generators ({active}) exceed platform bit width for indexing"
                )
            }
            Self::MismatchedGeneratorSet => {
                write!(
                    f,
                    "attempted to cast a CliffordNumber from an incompatible algebra"
                )
            }
        }
    }
}

#[cfg(feature = "std")]
impl Error for NumAnafisError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        None
    }
}

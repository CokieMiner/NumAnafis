//! Special-function namespace.
//!
//! The public special functions live as associated functions on [`SpecialFunc`],
//! implemented per topic in sibling files. A single import covers the whole API
//! with explicit call sites (`SpecialFunc::gamma`).

/// Zero-sized namespace for special mathematical functions.
///
/// Importing one item covers the public special-function API; call sites stay
/// explicit without polluting the local namespace.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct SpecialFunc;

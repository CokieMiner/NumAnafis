//! Integer arithmetic and representation types.

#![expect(
    clippy::trivially_copy_pass_by_ref,
    reason = "The arithmetic facade accepts borrowed operands independently of their representation"
)]

use core::cmp::Ordering;

// TODO(mp): replace fixed-width magnitude/root intermediates with backend
// square-root/cube-root remainder operations; retain borrowed arithmetic.
/// Native signed integer type.
pub type IntType = i64;

/// Zero-sized namespace for integer arithmetic over [`IntType`].
///
/// Addition, subtraction, multiplication, and negation wrap on overflow.
#[derive(Clone, Copy, Debug)]
#[non_exhaustive]
pub struct IntMath;

impl IntMath {
    /// Converts an `i64` into [`IntType`].
    #[inline]
    #[must_use]
    pub const fn from_i64(value: i64) -> IntType {
        value
    }

    /// Adds two integers, wrapping on overflow.
    #[inline]
    #[must_use]
    pub const fn add(lhs: &IntType, rhs: &IntType) -> IntType {
        (*lhs).wrapping_add(*rhs)
    }

    /// Subtracts two integers, wrapping on overflow.
    #[inline]
    #[must_use]
    pub const fn sub(lhs: &IntType, rhs: &IntType) -> IntType {
        (*lhs).wrapping_sub(*rhs)
    }

    /// Multiplies two integers, wrapping on overflow.
    #[inline]
    #[must_use]
    pub const fn mul(lhs: &IntType, rhs: &IntType) -> IntType {
        (*lhs).wrapping_mul(*rhs)
    }

    /// Wrapping negation: `i64::MIN` wraps to itself (modular semantics).
    #[inline]
    #[must_use]
    pub const fn wrapping_neg(value: &IntType) -> IntType {
        (*value).wrapping_neg()
    }
}

impl IntMath {
    /// Updates an integer accumulator, wrapping on overflow.
    pub const fn add_assign(lhs: &mut IntType, rhs: &IntType) {
        *lhs = Self::add(lhs, rhs);
    }
    /// Subtracts from an integer accumulator, wrapping on overflow.
    pub const fn sub_assign(lhs: &mut IntType, rhs: &IntType) {
        *lhs = Self::sub(lhs, rhs);
    }
    /// Exact length of an integer pair, when the result fits the integer backend.
    pub fn perfect_hypot(re: &IntType, im: &IntType) -> Option<IntType> {
        let squared = u128::from(re.unsigned_abs()).pow(2) + u128::from(im.unsigned_abs()).pow(2);
        let root = squared.isqrt();
        if root * root != squared {
            return None;
        }
        IntType::try_from(root).ok()
    }

    /// Compares two integers.
    #[inline]
    #[must_use]
    pub fn cmp(lhs: &IntType, rhs: &IntType) -> Ordering {
        lhs.cmp(rhs)
    }

    /// Returns the integer square root if `value` is a perfect square.
    #[inline]
    #[must_use]
    pub fn perfect_square(value: &IntType) -> Option<IntType> {
        {
            if *value < 0 {
                return None;
            }
            // The argument is nonnegative and its integer square root fits i64.
            let root = value.cast_unsigned().isqrt().cast_signed();
            (root * root == *value).then_some(root)
        }
    }

    /// Returns the integer cube root if `value` is a perfect cube.
    #[inline]
    #[must_use]
    pub fn perfect_cube(value: &IntType) -> Option<IntType> {
        Self::perfect_cube_unsigned(value.unsigned_abs()).map(|root| {
            // Cube roots of native integer magnitudes fit in the signed backend.
            let signed = root.cast_signed();
            if *value < 0 { -signed } else { signed }
        })
    }

    /// Exact root of an unsigned square denominator.
    #[must_use]
    pub const fn perfect_square_unsigned(value: u64) -> Option<u64> {
        let root = value.isqrt();
        if root * root == value {
            Some(root)
        } else {
            None
        }
    }
    /// Exact root of an unsigned cube denominator.
    #[must_use]
    pub fn perfect_cube_unsigned(value: u64) -> Option<u64> {
        if value == 0 {
            return Some(0);
        }
        let mut low = 0_u64;
        // A b-bit magnitude has cube root strictly below 2^ceil(b/3).
        let mut high = 1_u64 << value.bit_width().div_ceil(3);
        while low + 1 < high {
            let middle = low + ((high - low) >> 1);
            let cube = u128::from(middle) * u128::from(middle) * u128::from(middle);
            if cube <= u128::from(value) {
                low = middle;
            } else {
                high = middle;
            }
        }
        (u128::from(low) * u128::from(low) * u128::from(low) == u128::from(value)).then_some(low)
    }
}

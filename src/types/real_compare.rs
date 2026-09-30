//! Exact comparison and hashing of represented native real values.
//!
//! Finite values use a reduced odd fraction times a power of two. This permits
//! exact rational/float comparisons without allocating arbitrary integers.

use core::{
    cmp::Ordering,
    hash::{Hash, Hasher},
};

use super::{IntMath, RationalMath, Real};

// TODO(mp): generalize this canonical key to backend significands, exponents,
// and arbitrary integer magnitudes. Never compare exact values via f64.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct FiniteKey {
    negative: bool,
    numer: u64,
    denom: u64,
    shift: i32,
}

impl Real {
    /// Total order of represented real values, with NaNs ordered by their bits.
    /// Signed zeros compare equal; arithmetic comparisons keep NaNs unordered.
    #[must_use]
    pub fn total_cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .unwrap_or_else(|| match (self.is_nan(), other.is_nan()) {
                (true, true) => self.to_f64().to_bits().cmp(&other.to_f64().to_bits()),
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => Ordering::Equal,
            })
    }
}

impl PartialEq for Real {
    fn eq(&self, other: &Self) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

impl PartialOrd for Real {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.is_nan() || other.is_nan() {
            return None;
        }
        if self.is_infinite() || other.is_infinite() {
            return Some(infinity_class(self).cmp(&infinity_class(other)));
        }
        match (self, other) {
            (Self::Int(a), Self::Int(b)) => return Some(IntMath::cmp(a, b)),
            (Self::Rational(a), Self::Rational(b)) => return Some(RationalMath::cmp(a, b)),
            _ => {}
        }
        let a = finite_key(self);
        let b = finite_key(other);
        if a.negative != b.negative {
            return Some(if a.negative {
                Ordering::Less
            } else {
                Ordering::Greater
            });
        }
        let order = compare_positive(a, b);
        Some(if a.negative { order.reverse() } else { order })
    }
}

impl Hash for Real {
    fn hash<H: Hasher>(&self, state: &mut H) {
        if self.is_nan() {
            3_u8.hash(state);
            self.to_f64().to_bits().hash(state);
        } else if self.is_infinite() {
            infinity_class(self).hash(state);
        } else {
            1_u8.hash(state);
            finite_key(self).hash(state);
        }
    }
}

const fn infinity_class(value: &Real) -> u8 {
    if value.is_infinite() {
        if value.is_negative() { 0 } else { 2 }
    } else {
        1
    }
}

fn finite_key(component: &Real) -> FiniteKey {
    let (negative, numer, denom, shift) = match component {
        Real::Int(value) => (value.is_negative(), value.unsigned_abs(), 1, 0),
        Real::Rational(value) => (
            value.numer().is_negative(),
            value.numer().unsigned_abs(),
            *value.denom(),
            0,
        ),
        Real::Float(value) => {
            let bits = value.to_bits();
            let exponent = i32::try_from((bits >> 52) & 0x7ff).expect("eleven-bit exponent");
            let fraction = bits & ((1_u64 << 52) - 1);
            let (significand, shift) = if exponent == 0 {
                (fraction, -1074)
            } else {
                (fraction | (1_u64 << 52), exponent - 1075)
            };
            (bits >> 63 != 0, significand, 1, shift)
        }
    };
    if numer == 0 {
        return FiniteKey {
            negative: false,
            numer: 0,
            denom: 1,
            shift: 0,
        };
    }
    let nt = numer.trailing_zeros();
    let dt = denom.trailing_zeros();
    FiniteKey {
        negative,
        numer: numer >> nt,
        denom: denom >> dt,
        shift: shift + i32::try_from(nt).expect("at most 63 bits")
            - i32::try_from(dt).expect("at most 63 bits"),
    }
}

fn compare_positive(a: FiniteKey, b: FiniteKey) -> Ordering {
    let left = u128::from(a.numer) * u128::from(b.denom);
    let right = u128::from(b.numer) * u128::from(a.denom);
    if left == 0 || right == 0 {
        return left.cmp(&right);
    }
    let left_bits = i32::try_from(128 - left.leading_zeros()).expect("at most 128 bits");
    let right_bits = i32::try_from(128 - right.leading_zeros()).expect("at most 128 bits");
    match (left_bits + a.shift).cmp(&(right_bits + b.shift)) {
        Ordering::Equal => {
            if a.shift >= b.shift {
                (left
                    << u32::try_from(a.shift - b.shift).expect("equal lengths bound shift to 127"))
                .cmp(&right)
            } else {
                left.cmp(
                    &(right
                        << u32::try_from(b.shift - a.shift)
                            .expect("equal lengths bound shift to 127")),
                )
            }
        }
        order @ (Ordering::Less | Ordering::Greater) => order,
    }
}

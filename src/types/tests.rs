//! Numeric representation, arithmetic, and CAS inspection contracts.

extern crate std;

use core::{
    cmp::Ordering,
    hash::{Hash, Hasher},
};

use std::collections::hash_map::DefaultHasher;

use proptest::prelude::*;

use crate::{error::NumAnafisError, traits::Numeric};

use super::{Complex, IntMath, Number, RationalType, Real, c, i, n, r};

#[test]
fn integer_roots_and_lengths_use_exact_backend_magnitudes() {
    assert_eq!(IntMath::perfect_cube(&i64::MIN), Some(-(1_i64 << 21)));
    assert_eq!(IntMath::perfect_cube(&i64::MAX), None);
    assert_eq!(IntMath::perfect_cube_unsigned(u64::MAX), None);
    assert_eq!(n(i64::MIN).cbrt(), n(-(1_i64 << 21)));
    assert!(!n(i64::MIN).cbrt().is_float());
    assert_eq!(
        c(3_000_000_000_i64, 4_000_000_000_i64).abs(),
        n(5_000_000_000_i64)
    );
    assert!(!c(3_000_000_000_i64, 4_000_000_000_i64).abs().is_float());
    assert!(c(i64::MIN, i64::MIN).abs().is_finite());
}

#[test]
fn borrowed_real_accumulation_preserves_tiers_and_native_wrapping() {
    let mut integer = Real::from_int(i64::MAX);
    integer += &Real::one();
    assert_eq!(integer, Real::from_int(i64::MIN));
    integer -= &Real::one();
    assert_eq!(integer, Real::from_int(i64::MAX));
    let mut rational = Real::from_rational(RationalType::new(1, 2));
    rational += &Real::from_rational(RationalType::new(1, 2));
    assert!(matches!(rational, Real::Int(1)));
    rational -= &Real::from_float(0.5);
    assert!(matches!(rational, Real::Float(value) if value.to_bits() == 0.5_f64.to_bits()));
    let mut zero = Real::from_float(-0.0);
    zero += &Real::from_float(-0.0);
    assert!(zero.to_f64().is_sign_negative());
}

fn hash(value: &Number) -> u64 {
    let mut state = DefaultHasher::new();
    value.hash(&mut state);
    state.finish()
}

#[test]
fn represented_values_compare_without_float_rounding() {
    let a = n(1_i64 << 53);
    let b = n(a.to_f64());
    let c = n((1_i64 << 53) + 1);
    assert_eq!(a, b);
    assert_ne!(b, c);
    assert_eq!(a.total_cmp(&c), Ordering::Less);
    assert_ne!(r(1, 10), n(0.1));
    assert!(r(1, 10) < n(0.1));
    for (lhs, rhs) in [(n(1), n(1.0)), (n(0), n(-0.0)), (r(1, 2), n(0.5))] {
        assert_eq!(lhs, rhs);
        assert_eq!(hash(&lhs), hash(&rhs));
    }
    assert_ne!(Number::nan(), Number::nan());
    assert_eq!(Number::nan().partial_cmp(&n(1)), None);
}

#[test]
fn signed_unsigned_rationals_reduce_before_narrowing() {
    let a = RationalType::from_parts(1, 1_u64 << 32);
    let sum = &a + &a;
    assert_eq!(sum, RationalType::from_parts(1, 1_u64 << 31));
    assert_eq!(
        RationalType::new(1, i64::MIN),
        RationalType::from_parts(-1, 1_u64 << 63)
    );
    assert_eq!(
        RationalType::new(i64::MIN, i64::MIN),
        RationalType::from_integer(1)
    );
    assert_eq!(RationalType::from_parts(1, u64::MAX).to_integer(), 0);
    assert_eq!(r(3, 2).to_int(), 1);
    assert_eq!(r(-3, 2).to_int(), -1);
    assert_eq!(r(-3, 2).floor(), n(-2));
    assert_eq!(r(-3, 2).ceil(), n(-1));
    assert_eq!(r(-3, 2).round(), n(-2));
}

#[test]
fn exact_domain_promotion_and_cas_integrality() {
    assert_eq!(r(-4, 9).sqrt(), c(0, r(2, 3)));
    assert!(r(-4, 9).sqrt().im().is_rational());
    assert_eq!(r(2, 3).pow(&n(-3)), r(27, 8));
    assert_eq!(r(-8, 27).cbrt(), r(-2, 3));
    assert_eq!(i().pow(&n(4)), n(1));
    assert!(i().pow(&n(4)).is_int());
    assert!(n(2).is_int());
    assert!(n(2.0).is_int());
    assert!(!r(3, 2).is_int());
    assert!(!i().is_int());
    assert!(!n(f64::INFINITY).is_int());
    assert!(n(-2.0).pow(&n(0.5)).is_complex());
    assert_eq!(n(0).exp(), n(1));
    assert!(n(0).exp().is_int());
    assert!(n(0.0).exp().is_float());
    assert_eq!(n(1).pow(&n(i64::MIN)), n(1));
    assert!(n(0).pow(&n(-1)).is_infinite());
    assert!(n(i64::MIN).ln().is_finite());
    assert!(n(i64::MIN).sqrt().is_finite());
}

#[test]
fn integer_inspection_distinguishes_value_from_exactness_and_components() {
    assert!(n(2.0).is_int());
    assert!(n(2.0).is_float());
    assert!(!n(2.000_000_000_000_000_4).is_int());
    // The fraction has already rounded away; inspection cannot recover it.
    let rounded = n(9_007_199_254_740_992.0 + 0.5);
    assert!(rounded.is_int());
    assert!(rounded.is_float());
    for value in [Number::nan(), n(f64::INFINITY), n(f64::NEG_INFINITY)] {
        assert!(!value.is_int());
        assert!(!value.is_int_ring());
    }
    assert!(n(-0.0).is_int());
    assert!(n(3).is_int_ring());
    assert!(c(2, 3).is_int_ring());
    assert!(!c(2, 3).is_int());
    assert!(c(2.0, 3.0).is_int_ring());
    assert!(!c(2, r(3, 2)).is_int_ring());
    assert!(!r(3, 2).is_int_ring());
    let integer_fraction = RationalType::from_integer(5);
    assert_eq!(*integer_fraction.numer(), 5);
    assert_eq!(*integer_fraction.denom(), 1);
    assert_eq!(Number::from_rational(integer_fraction), n(5));
    let typed = n(5).extract_complex().expect("real is complex-compatible");
    assert_eq!(typed.re, Real::from_int(5));
    assert!(typed.im.is_zero());
}

#[test]
fn extrema_require_real_values_rather_than_an_artificial_complex_order() {
    assert_eq!(n(-3).max(&n(-2)), n(-2));
    assert_eq!(r(1, 2).min(&n(1)), r(1, 2));
    assert_eq!(n(2).clamp(&n(0), &n(1)), n(1));
    for value in [i(), c(2, 3), Number::nan()] {
        assert_eq!(
            value.partial_cmp(&value),
            (!value.is_nan()).then_some(Ordering::Equal)
        );
        assert!(value.max(&value).is_nan());
        assert!(value.min(&value).is_nan());
        assert!(value.max(&n(3)).is_nan());
        assert!(n(3).min(&value).is_nan());
        assert!(value.clamp(&n(0), &n(1)).is_nan());
    }
}

#[test]
fn complex_components_are_real_and_branch_zeros_survive() {
    let complex = Complex::new(n(2), n(3)).expect("real components");
    assert_eq!(complex.partial_cmp(&complex), Some(Ordering::Equal));
    assert_eq!(
        complex.partial_cmp(&Complex::new(n(2), n(4)).expect("real components")),
        None
    );
    assert_eq!(Complex::new(i(), n(0)), Err(NumAnafisError::ExpectedReal));
    let upper = c(-1, 0.0).sqrt();
    let lower = c(-1, -0.0).sqrt();
    assert_eq!(upper.im(), n(1));
    assert_eq!(lower.im(), n(-1));
    assert!(c(-1, -0.0).ln().im().is_negative());
    assert_eq!(c(1, i()), n(0));
    assert_eq!(c(3, 4).to_float().im(), n(4.0));
}

#[test]
fn scaled_complex_division_and_stable_near_zero_kernels() {
    let large = c(1e200, 1e200);
    assert_eq!(&large / &large, n(1));
    let tiny = c(1e-200, -1e-200);
    assert_eq!(&tiny / &tiny, n(1));
    let z = c(0.0, core::f64::consts::PI).expm1();
    assert!((z.re().to_f64() + 2.0).abs() < 1e-14);
    assert!((n(1e-20).expm1().to_f64() - 1e-20).abs() < 1e-35);
    assert!((n(1e-20).log1p().to_f64() - 1e-20).abs() < 1e-35);
    let tiny_complex = c(1e-20, 1e-20);
    let evaluated = tiny_complex.log1p();
    assert!((evaluated.re().to_f64() - 1e-20).abs() < 1e-35);
    assert!((evaluated.im().to_f64() - 1e-20).abs() < 1e-35);
}

proptest! {
    #[test]
    fn perfect_integer_cube_detection_matches_exact_arithmetic(root in -2_097_152_i64..2_097_152) {
        let cube = root.checked_pow(3).expect("bounded cube fits i64");
        prop_assert_eq!(IntMath::perfect_cube(&cube), Some(root));
        if root.unsigned_abs() > 1 {
            prop_assert_eq!(IntMath::perfect_cube(&(cube + 1)), None);
        }
    }
    #[test]
    fn borrowed_assignments_match_exact_arithmetic(a in -100_i64..100,b in -100_i64..100,d in 1_i64..100) {
        let first = Real::from_rational(RationalType::new(a, d));
        let second = Real::from_rational(RationalType::new(b, d));
        let expected_sum = &first + &second;
        let mut accumulator = first.clone();
        accumulator += &second;
        prop_assert_eq!(&accumulator, &expected_sum);
        accumulator -= &second;
        prop_assert_eq!(accumulator, first);
    }
    #[test]
    fn bounded_exact_arithmetic_obeys_ring_laws(a in -100_i64..100,b in -100_i64..100,c in -100_i64..100,d in 1_i64..100) {
        let first=r(a,d);let second=r(b,d);let third=r(c,d);
        prop_assert_eq!((&first+&second)+&third,&first+(&second+&third));
        prop_assert_eq!((&first*&second)*&third,&first*(&second*&third));
        prop_assert_eq!(&first*(&second+&third),&first*&second+&first*&third);
        prop_assert_eq!((&first+&second)-&second,first);
    }
    #[test]
    fn dyadic_values_share_equality_and_hash(numer in -100_000_i64..100_000,power in 0_u32..20) {
        let value=r(numer,1_i64<<power);
        let approximate=n(value.to_f64());
        prop_assert_eq!(&value,&approximate);
        prop_assert_eq!(hash(&value),hash(&approximate));
    }
    #[test]
    fn rational_order_and_projection_use_quotients(numer in -100_000_i64..100_000,denom in 1_i64..10_000) {
        let ratio=RationalType::new(numer,denom);
        prop_assert_eq!(ratio.to_integer(),numer.checked_div(denom).expect("bounded nonzero positive denominator"));
        let real=Real::from_rational(ratio);
        prop_assert!(real.floor()<=real);
        prop_assert!(real.ceil()>=real);
    }
}

#[test]
fn standard_conversions_preserve_values_and_reject_unsigned_capacity_loss() {
    let integer: Number = 42_i64.into();
    let approximate: Number = 42.0_f64.into();
    assert!(integer.is_int());
    assert!(approximate.is_float());
    assert_eq!(Number::try_from(42_u64).expect("native capacity"), integer);
    assert_eq!(
        Number::try_from(42_usize).expect("native capacity"),
        integer
    );
    assert!(
        Number::try_from(u64::MAX).is_err(),
        "conversion must not silently approximate"
    );
}

#[test]
fn primitive_literals_interoperate_without_extra_constructors() {
    assert_eq!(5 * i(), c(0, 5));
    assert_eq!(i() * 5, c(0, 5));
    assert_eq!(0.5 * i(), c(0, 0.5));
    assert_eq!(r(2, 3) * i(), c(0, r(2, 3)));
    assert_eq!(3 - r(1, 2), r(5, 2));
    assert_eq!(r(1, 2) - 3, r(-5, 2));
    assert_eq!(3 / &n(2), r(3, 2));
    assert_eq!(&n(2) / 3, r(2, 3));
}

//! Tests for 4-tier Number features.

use core::cmp::Ordering;
use num_anafis::{Numeric, c, i, n, r};

#[test]
fn test_float_tier_is_never_demoted() {
    let two = n(2.0);
    let zero = n(0.0);

    assert!(two.is_float());
    let root = two.sqrt();
    assert!(root.is_float());
    assert_eq!(root.to_f64(), 2.0_f64.sqrt());

    let exp_zero = zero.exp();
    assert!(exp_zero.is_float());
    assert_eq!(exp_zero.to_f64(), 1.0);

    let difference = &two - &two;
    assert!(difference.is_float());
    assert_eq!(difference.to_f64(), 0.0);

    let floored = n(2.7).floor();
    assert!(floored.is_float());
    assert_eq!(floored.to_f64(), 2.0);

    let powered = n(3.0).pow(&n(0));
    assert!(powered.is_float());
    assert_eq!(powered.to_f64(), 1.0);

    // Cancellation onto the real axis must not relabel a float as exact.
    let complex_cancellation = c(2.0, 1.0) - c(0.0, 1.0);
    assert!(!complex_cancellation.is_complex());
    assert!(complex_cancellation.is_float());
}

#[test]
fn test_exact_rational_roots_and_integer_powers_stay_exact() {
    let root = r(4, 9).sqrt();
    assert_eq!(root, r(2, 3));
    assert!(root.is_rational());

    let power = r(2, 3).pow(&n(-3));
    assert_eq!(power, r(27, 8));
    assert!(power.is_rational());

    let integer_power = n(3).pow(&n(4));
    assert_eq!(integer_power, n(81));
    assert!(integer_power.is_int());
}

#[test]
fn test_float_ieee_division_by_zero() {
    let positive_infinity = n(1.0) / n(0.0);
    let negative_infinity = n(-1.0) / n(0.0);
    let negative_zero_infinity = n(1.0) / n(-0.0);
    let indeterminate = n(0.0) / n(0.0);

    assert!(positive_infinity.is_float());
    assert_eq!(positive_infinity.to_f64(), f64::INFINITY);
    assert_eq!(negative_infinity.to_f64(), f64::NEG_INFINITY);
    assert_eq!(negative_zero_infinity.to_f64(), f64::NEG_INFINITY);
    assert!(indeterminate.to_f64().is_nan());
}

#[test]
fn test_float_nonfinite_operations_follow_ieee_results() {
    let inf_difference = n(f64::INFINITY) - n(f64::INFINITY);
    let zero_times_infinity = n(0.0) * n(f64::INFINITY);
    let exp_negative_infinity = n(f64::NEG_INFINITY).exp();
    let log_positive_zero = n(0.0).ln();

    assert!(inf_difference.to_f64().is_nan());
    assert!(zero_times_infinity.to_f64().is_nan());
    assert_eq!(exp_negative_infinity.to_f64(), 0.0);
    assert_eq!(log_positive_zero.to_f64(), f64::NEG_INFINITY);
}

#[test]
fn test_from_float_integer_boundary() {
    // Magnitudes at or above 2^63 are unrepresentable as i64 and must stay
    // Float (a saturating cast to Int(MAX) would silently corrupt the value).
    let two63 = f64::from_bits(0x43E0_0000_0000_0000); // Exactly 2^63; a decimal literal would trip lossy_float_literal.
    assert!(n(two63).is_float());
    assert!(n(1e19).is_float());
    assert_eq!(n(two63).to_f64(), two63);

    // Explicit Float construction stays approximate, even for integral values.
    let positive = n(42.0);
    let negative = n(-17.0);
    assert!(positive.is_float());
    assert!(negative.is_float());
    assert_eq!(positive.to_f64(), 42.0);
    assert_eq!(negative.to_f64(), -17.0);
    assert!((&positive + &n(1)).is_float());
    assert!(positive.floor().is_float());
    assert!(positive.to_float().is_float());
    assert!(n(42).is_int());
}

#[test]
fn test_constructors_and_promotion() {
    let int_val = n(42);
    assert!(int_val.is_int());

    let rat_val = r(22, 7);
    assert!(rat_val.is_rational());

    let flt_val = n(2.5);
    assert!(flt_val.is_float());

    let imag_val = 5 * i();
    assert!(imag_val.is_complex());
    assert_eq!(imag_val.re(), n(0));
    assert_eq!(imag_val.im(), n(5));

    let c_val = c(3, 4);
    assert!(c_val.is_complex());
    assert_eq!(c_val.re(), n(3));
    assert_eq!(c_val.im(), n(4));

    assert_eq!(i() * i(), n(-1));
}

#[test]
fn test_complex_branch_cuts() {
    assert_eq!(n(-1).sqrt(), i());

    let ln_neg1 = n(-1).ln();
    assert_eq!(ln_neg1.re(), n(0));
    assert!((ln_neg1.im().to_f64() - core::f64::consts::PI).abs() < 1e-12);
}

#[test]
fn test_complex_discretization_preserves_float_tier() {
    let z = c(3.7, 4.2);

    let z_floor = z.floor();
    assert_eq!(z_floor.re(), n(3));
    assert_eq!(z_floor.im(), n(4));
    assert!(z_floor.re().is_float());
    assert!(z_floor.im().is_float());

    let z_ceil = z.ceil();
    assert_eq!(z_ceil.re(), n(4));
    assert_eq!(z_ceil.im(), n(5));
    assert!(z_ceil.re().is_float());
    assert!(z_ceil.im().is_float());

    let z_round = z.round();
    assert_eq!(z_round.re(), n(4));
    assert_eq!(z_round.im(), n(4));

    let z_trunc = z.trunc();
    assert_eq!(z_trunc.re(), n(3));
    assert_eq!(z_trunc.im(), n(4));

    let z_im_zero = c(3.7, 0.4).floor();
    assert_eq!(z_im_zero, n(3));
    assert!(z_im_zero.is_float());
}

#[test]
fn test_sorting_order_is_separate_from_real_extrema() {
    let c1 = c(3, 4);
    let c2 = c(1, 2);

    assert!(c1.max(&c2).is_nan());
    assert!(c1.min(&c2).is_nan());
    assert_eq!(c1.total_cmp(&c2), Ordering::Greater);
    assert_eq!(c2.total_cmp(&c1), Ordering::Less);

    assert_eq!(n(7).max(&n(3)), n(7));
    assert_eq!(n(-5).max(&n(2)), n(2));
    assert_eq!(n(-5).min(&n(2)), n(-5));
}

#[test]
fn test_real_projection_to_f64() {
    let z = c(3, 4);
    assert_eq!(z.to_f64(), 3.0);
    assert_eq!(z.abs().to_f64(), 5.0);
}

#[test]
fn test_numeric_trait() {
    let z = c(0, core::f64::consts::PI / 2.0);
    let exp_z = Numeric::exp(&z);
    assert!((exp_z.re().to_f64()).abs() < 1e-12);
    assert!((exp_z.im().to_f64() - 1.0).abs() < 1e-12);
}

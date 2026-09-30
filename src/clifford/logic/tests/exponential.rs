//! Algebraic exponential shortcuts and their spectral fallback.

use alloc::vec;

use crate::{Numeric, types::Real};

use super::super::{Cga, CliffordNumber, GeneratorSet, SpectralDispatch, SpectralFn};

fn close(actual: &CliffordNumber, expected: &CliffordNumber) {
    assert!(actual.is_finite(), "exponential must be finite");
    for (a, b) in actual.coeffs_slice().iter().zip(expected.coeffs_slice()) {
        assert!(
            (a.to_f64() - b.to_f64()).abs() < 1e-12,
            "coefficient mismatch: {a} versus {b}"
        );
    }
}

#[test]
fn nilpotent_exponential_preserves_exact_coefficients() {
    let nilpotent = Cga::eps() * Real::from_int(3);
    let one = CliffordNumber::scalar(Cga::gens(), Real::one()).expect("real scalar");
    assert_eq!(nilpotent.exp(), &one + &nilpotent);
    assert_eq!(nilpotent.expm1(), nilpotent);
    assert!(
        nilpotent
            .exp()
            .coeffs_slice()
            .iter()
            .all(|value| matches!(value, Real::Int(_))),
        "nilpotent result must stay exact"
    );
    let shifted = one + nilpotent;
    let expected = shifted.clone() * Real::one().exp();
    close(&shifted.exp(), &expected);
}

#[test]
fn elliptic_and_hyperbolic_exponentials_use_the_signed_square() {
    for bivector in [Cga::qi(), Cga::sj()] {
        let input = bivector * Real::Float(0.3);
        let spectral = SpectralDispatch::apply_via_spectral(&input, &SpectralFn::Exp);
        close(&input.exp(), &spectral);
        let inverse = (-&input).exp();
        let one = CliffordNumber::scalar(Cga::gens(), Real::one()).expect("real scalar");
        close(&(&input.exp() * &inverse), &one);
    }
}

#[test]
fn nonsimple_bivector_retains_the_spectral_path() {
    let gens = GeneratorSet::new(vec![(1, 1), (2, 1), (3, 1), (4, 1)])
        .expect("valid generator labels, metrics, and count");
    let mut first = CliffordNumber::scalar(gens, Real::zero()).expect("scalar");
    let mut second = first.clone();
    *first.coeffs_mut_slice().get_mut(3).expect("e12 lane") = Real::Float(0.2);
    *second.coeffs_mut_slice().get_mut(12).expect("e34 lane") = Real::Float(0.3);
    let bivector = &first + &second;
    let square = &bivector * &bivector;
    assert!(
        !square.coeff(15).is_zero(),
        "B squared has a pseudoscalar part"
    );
    close(&bivector.exp(), &(&first.exp() * &second.exp()));
}

#[test]
fn expm1_shortcut_preserves_small_scalar_and_bivector_terms() {
    let one = CliffordNumber::scalar(Cga::gens(), Real::Float(1e-20)).expect("real scalar");
    let input = one + Cga::qi() * Real::Float(1e-20);
    let result = input.expm1();
    assert!(
        (result.coeff(0).to_f64() - 1e-20).abs() < 1e-35,
        "tiny scalar expm1 term must survive cancellation"
    );
    assert!(
        (result.coeff(6).to_f64() - input.coeff(6).to_f64()).abs() < 1e-35,
        "tiny bivector term must survive cancellation"
    );
}

#[test]
fn single_blade_exponential_avoids_unnecessary_squaring() {
    let coefficient = Real::from_rational(crate::RationalType::from_parts(1, u64::MAX));
    let input = Cga::qi() * coefficient.clone();
    let result = input.exp();
    assert!(
        result.is_finite(),
        "large rational denominator must not be squared"
    );
    assert!(
        (result.coeff(6).to_f64() + coefficient.to_f64()).abs() < 1e-34,
        "tiny rotation coefficient must survive"
    );
    let large = Cga::qi() * Real::Float(1e200);
    assert!(
        large.exp().is_finite(),
        "finite rotation magnitude must not overflow when squared"
    );
}

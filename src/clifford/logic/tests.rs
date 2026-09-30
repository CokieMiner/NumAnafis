//! Tests for Clifford constructors and algebraic identities.

#[path = "tests/constructors.rs"]
mod constructors;

#[path = "tests/exponential.rs"]
mod exponential;

#[path = "tests/api.rs"]
mod api;

use super::{
    Cga, CliffordNumber, Cmplx, EigenSolver, GeneratorSet, Jet, LargeEmbed, Mat2C, MatC,
    SmallEmbed, SpectralFn,
};
use crate::{
    NumAnafisError,
    traits::Numeric,
    types::{Number, RationalType, Real, r},
};

use alloc::vec::Vec;

fn test_scalar(val: Real) -> CliffordNumber {
    CliffordNumber::scalar(Cga::gens(), val).expect("scalar creation is infallible")
}

#[test]
fn coefficient_distance_preserves_exact_tolerance_boundaries() {
    let zero = CliffordNumber::zero(Cga::gens()).expect("valid algebra");
    let approximate = test_scalar(Real::from_float(0.1));
    let exact_bound = CliffordNumber::scalar(Cga::gens(), r(1, 10)).expect("real coefficient");
    assert!(!approximate.approx_eq_number(&zero, &exact_bound));
    assert!(approximate.approx_eq_number(&zero, &approximate));
    assert!(!test_scalar(Real::nan()).approx_eq_number(&zero, &exact_bound));
}

#[test]
fn test_cga_base_squares() {
    // e1² = e2² = e3² = e+² = +1
    assert_eq!(&Cga::e1() * &Cga::e1(), test_scalar(Real::one()));
    assert_eq!(&Cga::e2() * &Cga::e2(), test_scalar(Real::one()));
    assert_eq!(&Cga::e3() * &Cga::e3(), test_scalar(Real::one()));
    assert_eq!(&Cga::e_plus() * &Cga::e_plus(), test_scalar(Real::one()));

    // e−² = -1
    assert_eq!(&Cga::e_minus() * &Cga::e_minus(), test_scalar(-Real::one()));
}

#[test]
fn test_derived_identities() {
    // ε² = 0 (nilpotent, like dual numbers)
    assert!(
        (&Cga::eps() * &Cga::eps())
            .coeffs_slice()
            .iter()
            .all(Real::is_zero)
    );

    // ci² = -1 (complex imaginary, like complex numbers)
    let ci_sq = &Cga::ci() * &Cga::ci();
    assert_eq!(ci_sq.coeffs_slice()[0], -Real::one());

    // sj² = +1 (split-complex unit)
    let sj_sq = &Cga::sj() * &Cga::sj();
    assert_eq!(sj_sq.coeffs_slice()[0], Real::one());

    // Quaternions: i² = j² = k² = ijk = -1
    assert_eq!((&Cga::qi() * &Cga::qi()).coeffs_slice()[0], -Real::one());
    assert_eq!((&Cga::qj() * &Cga::qj()).coeffs_slice()[0], -Real::one());
    assert_eq!((&Cga::qk() * &Cga::qk()).coeffs_slice()[0], -Real::one());

    let ijk = &(&Cga::qi() * &Cga::qj()) * &Cga::qk();
    assert_eq!(ijk.coeffs_slice()[0], -Real::one());
}

#[test]
fn test_cga_central_complex_unit() {
    let imaginary = Cga::ci();
    assert_eq!(imaginary, Cga::pseudo5d());
    assert_eq!(&imaginary * &imaginary, test_scalar(-Real::one()));
    for generator in [
        Cga::e1(),
        Cga::e2(),
        Cga::e3(),
        Cga::e_plus(),
        Cga::e_minus(),
    ] {
        assert_eq!(&imaginary * &generator, &generator * &imaginary);
    }

    // The Euclidean pseudoscalar is a different complex structure: it
    // anticommutes with the two additional conformal generators.
    let euclidean_imaginary = Cga::pseudo3d();
    for generator in [Cga::e_plus(), Cga::e_minus()] {
        let left = &euclidean_imaginary * &generator;
        let right = &generator * &euclidean_imaginary;
        assert_eq!(left, -right);
        assert!(left.coeffs_slice().iter().any(|value| !value.is_zero()));
    }

    // The central embedding preserves ordinary complex multiplication:
    // (2 + 3i)(4 - 5i) = 23 + 2i.
    let a = test_scalar(Real::from_int(2)) + &imaginary * Real::from_int(3);
    let b = test_scalar(Real::from_int(4)) - &imaginary * Real::from_int(5);
    let expected = test_scalar(Real::from_int(23)) + &imaginary * Real::from_int(2);
    assert_eq!(a * b, expected);
}

#[test]
fn test_conformal_points() {
    // eₒ² = 0, e∞² = 0
    assert!(
        (&Cga::orig() * &Cga::orig())
            .coeffs_slice()
            .iter()
            .all(Real::is_zero)
    );
    assert!(
        (&Cga::inf() * &Cga::inf())
            .coeffs_slice()
            .iter()
            .all(Real::is_zero)
    );

    // eₒ · e∞ + e∞ · eₒ = -2 (since eₒ · e∞ = -1)
    let op = &Cga::orig() * &Cga::inf();
    let po = &Cga::inf() * &Cga::orig();
    let dot = &op + &po;
    assert_eq!(dot.coeffs_slice()[0], -(Real::one() + Real::one()));
}

#[test]
fn test_spectral_round_trip() {
    let val_2 = test_scalar(Real::one() + Real::one());
    let val_1 = test_scalar(Real::one());
    let val_05 = test_scalar(Real::one() / (Real::one() + Real::one()));

    let mv = &val_2 + &(&Cga::e1() * &val_1) + &(&Cga::e2() * &val_05);

    // Compute exp(mv) and then ln(exp(mv)) to round-trip
    let exp_mv = mv.exp();
    let round_trip = exp_mv.ln();

    let eps = Real::epsilon();
    let mut multiplier = Real::one();
    for _ in 0..49 {
        multiplier = &multiplier + &Real::one();
    }
    let tol = &eps * &multiplier;

    let limit = mv.blade_count();
    for i in 0..limit {
        let diff = (&round_trip.coeffs_slice()[i] - &mv.coeffs_slice()[i]).abs();
        assert!(
            diff.total_cmp(&tol) == core::cmp::Ordering::Less,
            "Real at {} differs: got {}, expected {}, diff: {}, tol: {}",
            i,
            round_trip.coeffs_slice()[i],
            mv.coeffs_slice()[i],
            diff,
            tol
        );
    }
}

#[test]
fn test_scalar_mv_add_sub() {
    let one = Real::one();

    // MV + scalar: only the scalar blade moves.
    let sum_r = &Cga::e1() + &one;
    assert_eq!(sum_r.coeffs_slice()[0], one);
    assert_eq!(sum_r.coeffs_slice()[1], one);
    assert!(sum_r.coeffs_slice()[2].is_zero());

    // Real + MV agrees (addition commutes for grade-0).
    let sum_l = &one + &Cga::e1();
    assert_eq!(sum_l, sum_r);

    // MV - scalar: only the scalar blade moves.
    let diff = &Cga::e1() - &one;
    assert_eq!(diff.coeffs_slice()[0], -one.clone());
    assert_eq!(diff.coeffs_slice()[1], one);
}

#[test]
fn test_translator_inverse() {
    // T = 1 - e1*e_inf is a CGA translator: N = e1*e_inf is nilpotent
    // (N^2 = 0 since e_inf^2 = 0), so T^-1 = 1 + N exactly. Inversion
    // goes through LU (no eigenvectors), which handles the defective
    // matrix image that eigendecomposition cannot.
    let n = &Cga::e1() * &Cga::inf();
    assert!((&n * &n).coeffs_slice().iter().all(Real::is_zero));
    let t = &test_scalar(Real::one()) - &n;
    let inv = t.geometric_inverse();
    let expect = &test_scalar(Real::one()) + &n;

    let tol = Real::epsilon() * Real::from_int(1024);
    let limit = t.blade_count();
    for i in 0..limit {
        let diff = (&inv.coeffs_slice()[i] - &expect.coeffs_slice()[i]).abs();
        assert!(
            diff.total_cmp(&tol) == core::cmp::Ordering::Less,
            "Blade {} inverts poorly: got {}, expected {}, diff: {}",
            i,
            inv.coeffs_slice()[i],
            expect.coeffs_slice()[i],
            diff
        );
    }
}

// Representation contracts.

#[test]
fn test_blade_count_is_power_of_two() {
    for (n, ids) in [
        (0, alloc::vec![]),
        (1, alloc::vec![(300, 1)]),
        (3, alloc::vec![(300, 1), (301, 1), (302, -1)]),
        (
            5,
            alloc::vec![(300, 1), (301, 1), (302, 1), (303, 1), (304, -1)],
        ),
    ] {
        let gens = GeneratorSet::new(ids).expect("valid generator labels, metrics, and count");
        let mv = CliffordNumber::zero_unchecked(gens);
        assert_eq!(mv.blade_count(), 1 << n);
    }
}

#[test]
#[should_panic(expected = "out of range")]
fn test_blade_access_rejects_padding() {
    let gens = GeneratorSet::new(alloc::vec![(300, 1)])
        .expect("valid generator labels, metrics, and count");
    let mv = CliffordNumber::zero_unchecked(gens);
    assert!(mv.coeff(2).is_zero());
}

// Small-embedding fidelity contracts.

fn roundtrip_small(metrics: &[i8]) -> Option<Real> {
    let entries: Vec<(u32, i8)> = metrics
        .iter()
        .enumerate()
        .map(|(i, m)| {
            (
                u32::try_from(400).expect("test id fits in u32")
                    + u32::try_from(i).expect("test index fits in u32"),
                *m,
            )
        })
        .collect();
    let gens = GeneratorSet::new(entries).expect("valid generator labels, metrics, and count");
    let basis = SmallEmbed::compute_blade_basis(&gens)?;
    let n = gens.len();
    let mut worst = Real::zero();
    for blade in 0..(1 << n) {
        let mut mv = CliffordNumber::zero_unchecked(gens.clone());
        mv.coeffs_mut_slice()[blade] = Real::one();
        let rt = SmallEmbed::from_matrix_with_gens(
            &SmallEmbed::to_matrix_with_basis(&mv, &basis),
            &gens,
            &basis,
        );
        for b in 0..(1 << n) {
            let want = if b == blade {
                Real::one()
            } else {
                Real::zero()
            };
            let d = (rt.coeff(b) - &want).abs();
            if d.total_cmp(&worst) == core::cmp::Ordering::Greater {
                worst = d;
            }
        }
    }
    Some(worst)
}

#[test]
fn test_small_embedding_roundtrip() {
    // Faithful signatures round-trip exactly; unfaithful ones decline
    // (None) so callers fall back to the large path.
    for metrics in [
        &[1][..],
        &[-1][..],
        &[0][..],
        &[1, 1][..],
        &[1, -1][..],
        &[-1, -1][..],
        &[0, -1][..],
        &[-1, 0][..],
        &[1, 1, 1][..],
        &[1, -1, -1][..],
        &[-1, -1, 1][..],
        &[-1, 1, -1][..],
    ] {
        let err = roundtrip_small(metrics).expect("faithful signature must embed");
        assert!(err.is_zero(), "{metrics:?} leaks: {err}");
    }
    for metrics in [
        &[0, 1][..],
        &[1, 0][..],
        &[-1, -1, -1][..],
        &[-1, 1, 1][..],
        &[1, -1, 1][..],
        &[1, 1, -1][..],
    ] {
        assert!(
            roundtrip_small(metrics).is_none(),
            "{metrics:?} must decline (unfaithful)"
        );
    }
}

// Eigensolver contracts.

fn mat2c_recon_error(m: &Mat2C) -> Option<Cmplx> {
    let (u, uinv, l1, l2) = SmallEmbed::eigendecompose(m)?;
    let recon = u
        .mul(&Mat2C::new(l1, Cmplx::zero(), Cmplx::zero(), l2))
        .mul(&uinv);
    let mut worst = Cmplx::zero();
    for (x, y) in [
        (&recon.a, &m.a),
        (&recon.b, &m.b),
        (&recon.c, &m.c),
        (&recon.d, &m.d),
    ] {
        let d = x.sub(y).abs();
        if d.total_cmp(&worst.abs()) == core::cmp::Ordering::Greater {
            worst = x.sub(y);
        }
    }
    Some(worst)
}

#[test]
fn test_2x2_diagonal_order() {
    // eigendecompose(diag(1,2)) must pair U = I with eigenvalues in place.
    let d = Mat2C::new(
        Cmplx::new(Real::one(), Real::zero()),
        Cmplx::new(Real::zero(), Real::zero()),
        Cmplx::new(Real::zero(), Real::zero()),
        Cmplx::new(Real::one() + Real::one(), Real::zero()),
    );
    let err = mat2c_recon_error(&d).expect("diagonal always diagonalizes");
    assert!(err.abs().is_zero(), "diagonal order wrong: {err:?}");
}

#[test]
fn test_2x2_lower_triangular() {
    // [[1,0],[1,2]] has distinct eigenvalues; the λ=d eigenvector (0,1)
    // must not be rejected.
    let lt = Mat2C::new(
        Cmplx::new(Real::one(), Real::zero()),
        Cmplx::new(Real::zero(), Real::zero()),
        Cmplx::new(Real::one(), Real::zero()),
        Cmplx::new(Real::one() + Real::one(), Real::zero()),
    );
    let err = mat2c_recon_error(&lt).expect("distinct eigenvalues diagonalize");
    assert!(err.abs().is_zero(), "lower-triangular wrong: {err:?}");
}

// Real-domain contracts.

#[test]
fn test_scalar_complex_continuation() {
    let gens = Cga::gens();
    let neg1 = CliffordNumber::scalar(gens.clone(), -Real::one()).expect("scalar");
    let root = neg1.sqrt();
    assert_eq!(root, Cga::ci());
    assert_eq!(root.generator_set(), &gens);
    let logarithm = neg1.ln();
    assert_eq!(logarithm.coeff(0), &Real::zero());
    assert_eq!(logarithm.coeff(31), &Real::pi());
    let imaginary = CliffordNumber::scalar(gens.clone(), Number::i()).expect("central scalar");
    let expected = Number::i().exp();
    assert_eq!(
        Number::from(imaginary.exp())
            .extract_complex()
            .expect("complex scalar"),
        expected.extract_complex().expect("complex value")
    );
    let input = CliffordNumber::scalar(
        gens,
        Real::from_rational(crate::types::RationalType::new(-4, 9)),
    )
    .expect("scalar");
    assert_eq!(
        input.sqrt(),
        Cga::ci() * Real::from_rational(RationalType::from_parts(2, 3))
    );
    let unsupported = CliffordNumber::scalar(GeneratorSet::empty(), -Real::one())
        .expect("real scalar")
        .sqrt();
    assert!(!unsupported.is_finite());
}

// Large-pipeline contracts.

#[test]
fn test_large_embedding_roundtrip() {
    // Every signature round-trips through the large path, including ones
    // the small path declines.
    for metrics in [
        &[0, 1][..],
        &[1, 0][..],
        &[-1, -1, -1][..],
        &[0, 0][..],
        &[1, 1, 1][..],
        &[1, 1, 1, 1, 1][..],
        &[1, 1, 1, 1, -1][..],
    ] {
        let entries: Vec<(u32, i8)> = metrics
            .iter()
            .enumerate()
            .map(|(i, m)| {
                (
                    u32::try_from(500).expect("test id fits in u32")
                        + u32::try_from(i).expect("test index fits in u32"),
                    *m,
                )
            })
            .collect();
        let gens = GeneratorSet::new(entries).expect("valid generator labels, metrics, and count");
        let basis = LargeEmbed::blade_basis(&gens).expect("large basis exists");
        let n = gens.len();
        for blade in 0..(1 << n) {
            let mut mv = CliffordNumber::zero_unchecked(gens.clone());
            mv.coeffs_mut_slice()[blade] = Real::one();
            let rt = LargeEmbed::from_matrix(&LargeEmbed::to_matrix(&mv, &basis), &gens, &basis);
            for b in 0..(1 << n) {
                let want = if b == blade {
                    Real::one()
                } else {
                    Real::zero()
                };
                assert_eq!(rt.coeff(b), &want, "{metrics:?} blade {blade} leaks");
            }
        }
    }
}

#[test]
fn test_qr_reconstruction() {
    // Real symmetric cases reconstruct to float precision; eigenvalues exact.
    let mut m = MatC::zeros(4);
    for i in 0..4 {
        m.set(i, i, Cmplx::new(Real::one() + Real::one(), Real::zero()));
        if i > 0 {
            m.set(i, i - 1, Cmplx::new(Real::one(), Real::zero()));
            m.set(i - 1, i, Cmplx::new(Real::one(), Real::zero()));
        }
    }
    let (q, t) = EigenSolver::schur_decomposition(&m).expect("schur converges");
    let recon = q.mul(&t).mul(&q.dagger());
    let tol = Real::epsilon() * Real::from_int(10000);
    for i in 0..4 {
        for j in 0..4 {
            let d = m.get(i, j).sub(recon.get(i, j)).abs();
            assert!(
                d.total_cmp(&tol) == core::cmp::Ordering::Less,
                "QR residual at ({i},{j}): {d:?}"
            );
        }
    }
}

#[test]
fn test_sinc_jet_composes() {
    // sinc(2ε) = 1 − (2/3)ε² + …, not the bare series 1 − (1/6)ε².
    let two = Real::one() + Real::one();
    let x = Jet::new(alloc::vec![
        Cmplx::new(Real::zero(), Real::zero()),
        Cmplx::new(two, Real::zero()),
        Cmplx::new(Real::zero(), Real::zero()),
    ]);
    let scalar_ratio = SpectralFn::Sinc.eval_jet(&x);
    assert!(scalar_ratio.coeffs[0].0.is_one());
    assert!(scalar_ratio.coeffs[1].is_zero());
    let want = Real::from_int(-2) / Real::from_int(3);
    let tol = Real::epsilon() * Real::from_int(1000);
    let diff = (scalar_ratio.coeffs[2].0.clone() - want).abs();
    assert!(
        diff.total_cmp(&tol) == core::cmp::Ordering::Less,
        "sinc composition wrong: {diff:?}"
    );
}

#[test]
fn test_cbrt_branches_agree() {
    // Real and Jet paths take the real branch on real input: cbrt(-8) = -2.
    assert_eq!(Real::from_int(-8).cbrt(), Real::from_int(-2));
    let x = Jet::variable(&Cmplx::new(Real::from_int(-8), Real::zero()), 2);
    let scalar_ratio = SpectralFn::Cbrt.eval_jet(&x);
    assert_eq!(scalar_ratio.coeffs[0].0, Real::from_int(-2));
    assert!(scalar_ratio.coeffs[0].1.is_zero());
}

#[test]
fn test_sqrt_outside_algebra_is_nan() {
    // √e ∉ Cl(1,0) (no real a + be squares to e): NaN, not (1+e)/2.
    let gens = GeneratorSet::new(alloc::vec![(600, 1)])
        .expect("valid generator labels, metrics, and count");
    let e = CliffordNumber::generator(gens, 0).expect("generator");
    let scalar = e.sqrt();
    assert!(scalar.coeff(0).is_nan());
}

#[test]
fn test_expm1_tiny() {
    // expm1(1e-20) ≈ 1e-20, not exactly 0 (cancellation-free kernel).
    let tiny = Real::from_int(1) / Real::from_int(10).pow(&Real::from_int(20));
    let close = (tiny.expm1() - tiny).abs();
    let tol = Real::epsilon() * Real::from_int(100);
    assert!(
        close.total_cmp(&tol) == core::cmp::Ordering::Less,
        "expm1 cancellation: {close:?}"
    );
}

#[test]
fn central_matrix_units_and_function_outputs_agree() {
    for metrics in [
        alloc::vec![-1],
        alloc::vec![1, 1, 1],
        alloc::vec![1, 1, 1, 1, -1],
    ] {
        let gens = GeneratorSet::new(
            metrics
                .iter()
                .enumerate()
                .map(|(index, &metric)| (u32::try_from(index).expect("small index"), metric))
                .collect(),
        )
        .expect("valid generator labels, metrics, and count");
        let blade = gens.complex_blade().expect("central unit");
        let unit = CliffordNumber::scalar(gens.clone(), Number::i()).expect("central scalar");
        let basis = LargeEmbed::blade_basis(&gens).expect("faithful embedding");
        let mapped = LargeEmbed::to_matrix(&unit, &basis);
        assert_eq!(mapped.dim, 1 << ((gens.len() - 1) >> 1));
        for row in 0..mapped.dim {
            for column in 0..mapped.dim {
                assert_eq!(
                    mapped.get(row, column),
                    &if row == column {
                        Cmplx::i()
                    } else {
                        Cmplx::zero()
                    }
                );
            }
        }
        assert_eq!(Number::from(unit.clone()), Number::i());
        assert_eq!(Number::from(unit.clone()).conj(), -Number::i());
        assert_eq!(Number::from(unit.clone()).norm_sq(), Number::one());
        assert_eq!(unit.coeff(blade), &Real::one());
        let vector = CliffordNumber::generator(gens.clone(), 0).expect("generator");
        let shifted = Number::from(vector.clone()) + Number::i();
        let expected = CliffordNumber::scalar(gens, Number::i()).expect("central") + vector;
        assert_eq!(shifted, Number::from(expected));
    }
    let gens = Cga::gens();
    // e1 has eigenvalues +/-1. sqrt maps the negative eigenspace into I5.
    let root = Cga::e1().sqrt();
    let expected = (test_scalar(Real::one()) + Cga::e1())
        * Real::from_rational(RationalType::from_parts(1, 2))
        + Cga::ci()
            * (test_scalar(Real::one()) - Cga::e1())
            * Real::from_rational(RationalType::from_parts(1, 2));
    assert!(root.approx_eq_number(&expected, &test_scalar(Real::from_float(1e-12))));
    assert!(root.coeffs_slice().iter().all(Real::is_finite));
    assert_eq!(
        CliffordNumber::scalar(GeneratorSet::empty(), Number::i()),
        Err(NumAnafisError::UnsupportedComplexStructure)
    );
    assert_eq!(gens.complex_blade(), Some(31));
}

#[test]
fn scalar_meaning_survives_different_clifford_contexts() {
    let gens = GeneratorSet::new(alloc::vec![(0, 1), (1, 1), (2, 1)])
        .expect("valid generator labels, metrics, and count");
    let i3: Number = CliffordNumber::scalar(gens, Number::i())
        .expect("complex structure")
        .into();
    let i5: Number = Cga::ci().into();
    assert_eq!(i3, i5);
    assert_eq!(&i3 + &i5, Number::from(Cga::ci() * 2_i64));
    assert_eq!(&i3 * &i5, Number::neg_one());
    let vector: Number = Cga::e1().into();
    assert_eq!(&vector * &i3, &vector * &Number::i());
    let negative = test_scalar(Real::from_int(-4));
    let half = test_scalar(Real::from_rational(RationalType::from_parts(1, 2)));
    assert_eq!(negative.pow(&half), Cga::ci() * 2_i64);
    let two = test_scalar(Real::from_int(2));
    assert_eq!(Cga::qi().pow(&two), test_scalar(Real::neg_one()));
    assert!(Cga::qi().pow(&two).coeff(0).is_int());
}

#[test]
fn lower_jordan_exponential_is_identity_plus_nilpotent() {
    let gens = GeneratorSet::new(alloc::vec![(0, 1), (1, -1)])
        .expect("valid generator labels, metrics, and count");
    let a = CliffordNumber::generator(gens.clone(), 0).expect("generator");
    let b = CliffordNumber::generator(gens.clone(), 1).expect("generator");
    let nilpotent = (a - b) * Real::from_rational(RationalType::from_parts(1, 2));
    assert!(nilpotent.geometric_mul(&nilpotent).is_zero());
    let expected = CliffordNumber::scalar(gens, Real::one()).expect("identity") + &nilpotent;
    assert!(
        nilpotent.exp().approx_eq_number(
            &expected,
            &CliffordNumber::scalar(expected.generator_set().clone(), Real::from_float(1e-12))
                .expect("tolerance")
        )
    );
}

#[test]
fn jet_products_stay_at_requested_degree_and_preserve_low_orders() {
    let variable = Jet::variable(&Cmplx::zero(), 4);
    let cube = variable.mul(&variable).mul(&variable);
    assert_eq!(cube.degree(), 4);
    assert_eq!(cube.coeffs[3], Cmplx::one());
    assert!(
        cube.coeffs
            .iter()
            .enumerate()
            .all(|(index, value)| index == 3 || value.is_zero())
    );
    let series = SpectralFn::Sinc.eval_jet(&Jet::variable(&Cmplx::zero(), 64));
    assert_eq!(series.degree(), 64);
    assert!((series.coeffs[22].0.to_f64() + 1.0 / libm::tgamma(24.0)).abs() < 1e-35);
    let tiny = Jet::variable(&Cmplx::new(Real::from_float(1e-20), Real::zero()), 2);
    assert!((SpectralFn::Expm1.eval_jet(&tiny).coeffs[0].0.to_f64() - 1e-20).abs() < 1e-35);
    assert!((SpectralFn::Log1p.eval_jet(&tiny).coeffs[0].0.to_f64() - 1e-20).abs() < 1e-35);
}

#[test]
fn fixed_signature_products_match_dynamic_blade_algebra() {
    use super::FastClifford;
    fn check<const P: usize, const Q: usize, const R: usize>() {
        let count = 1 << (P + Q + R);
        for a in 0..count {
            for b in 0..count {
                let mut lhs = FastClifford::<P, Q, R>::zero();
                lhs.coeffs[a] = Real::one();
                let mut rhs = FastClifford::<P, Q, R>::zero();
                rhs.coeffs[b] = Real::one();
                let expected = CliffordNumber::from(&lhs) * CliffordNumber::from(&rhs);
                let mut assigned = lhs.clone();
                assigned *= &rhs;
                assert_eq!(CliffordNumber::from(&assigned), expected);
                let mut owned_assigned = lhs.clone();
                owned_assigned *= rhs.clone();
                assert_eq!(CliffordNumber::from(&owned_assigned), expected);
                assert_eq!(CliffordNumber::from(&(lhs * rhs)), expected);
            }
        }
    }
    check::<4, 1, 0>();
    check::<1, 1, 1>();
    check::<0, 1, 0>();
    check::<0, 0, 0>();
}

#[test]
fn transitive_cluster_ids_are_carried_into_parlett_blocks() {
    use super::large_spectral::{apply_ast_schur, reorder_schur};
    let mut t = MatC::zeros(4);
    for (index, value) in [0.0, 2e-7, 1e-7, 10.0].iter().enumerate() {
        t.set(
            index,
            index,
            Cmplx::new(Real::from_float(*value), Real::zero()),
        );
    }
    let blocks = reorder_schur(&mut t, &mut MatC::identity(4));
    assert_eq!(blocks, alloc::vec![(0, 2), (3, 3)]);
    let evaluated =
        apply_ast_schur(&t, &SpectralFn::Exp, &blocks).expect("cluster series converges");
    for index in 0..4 {
        assert!(
            (evaluated.get(index, index).0.to_f64() - libm::exp(t.get(index, index).0.to_f64()))
                .abs()
                < 1e-10
        );
    }
    // A finite Taylor budget must not silently return a divergent log series.
    let mut spread = MatC::zeros(2);
    spread.set(0, 0, Cmplx::one());
    spread.set(1, 1, Cmplx::new(Real::from_float(3.0), Real::zero()));
    assert!(apply_ast_schur(&spread, &SpectralFn::Ln, &[(0, 1)]).is_none());
}

#[test]
fn number_conjugation_is_a_uniform_complex_anti_involution() {
    for gens in [
        GeneratorSet::new(alloc::vec![(0, 1), (1, 1), (2, 1)])
            .expect("valid generator labels, metrics, and count"),
        Cga::gens(),
    ] {
        let first: Number = CliffordNumber::generator(gens.clone(), 0)
            .expect("generator")
            .into();
        let second: Number = CliffordNumber::generator(gens, 1)
            .expect("generator")
            .into();
        let a = &first * Number::from_int(2) + Number::i();
        let b = &second + Number::i() * Number::from_int(3);
        assert_eq!((&a + &b).conj(), a.conj() + b.conj());
        assert_eq!((&a * &b).conj(), b.conj() * a.conj());
        assert_eq!(a.conj().conj(), a);
        let expected = if first.as_clifford().expect("vector").n_generators() == 3 {
            5
        } else {
            -3
        };
        assert_eq!(a.norm_sq(), Number::from_int(expected));
    }
}

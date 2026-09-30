//! CGA catalog contracts and independent multiplication identities.

use alloc::{vec, vec::Vec};

use proptest::prelude::*;

use crate::{
    clifford::{Cga, CliffordNumber},
    types::{Number, Real},
};

#[test]
fn oriented_volume_and_quaternion_definitions() {
    let (x, y, z) = (Cga::e1(), Cga::e2(), Cga::e3());
    let i3 = &x * &y * &z;
    assert_eq!(Cga::pseudo3d(), i3);
    assert_eq!(Cga::pseudo3d().grade(3), Cga::pseudo3d());
    assert_eq!(Cga::pseudo5d(), &i3 * &Cga::e_plus() * &Cga::e_minus());
    assert_eq!(Cga::pseudo5d().grade(5), Cga::pseudo5d());
    assert_eq!(Cga::qi(), -(&y * &z));
    assert_eq!(Cga::qj(), -(&z * &x));
    assert_eq!(Cga::qk(), -(&x * &y));
    for (unit, axis) in [(Cga::qi(), x), (Cga::qj(), y), (Cga::qk(), z)] {
        assert_eq!(unit, -(&Cga::pseudo3d() * &axis));
        assert_eq!(unit.grade(2), unit);
    }
}

#[test]
fn hamilton_multiplication_and_conjugation() {
    let (i, j, k) = (Cga::qi(), Cga::qj(), Cga::qk());
    for (a, b, expected) in [(&i, &j, &k), (&j, &k, &i), (&k, &i, &j)] {
        assert_eq!(a * b, *expected);
        assert_eq!(b * a, -expected);
    }
    let q = scalar(2) + &i * 3_i64 + &j * 4_i64 + &k * 5_i64;
    let conjugate = scalar(2) - &i * 3_i64 - &j * 4_i64 - &k * 5_i64;
    assert_eq!(q.reverse(), conjugate);
    assert_eq!(q.clifford_conjugate(), conjugate);
    assert_eq!(&q * &conjugate, scalar(54));
}

#[test]
fn split_complex_and_dual_arithmetic() {
    let j = Cga::sj();
    let eps = Cga::eps();
    assert_eq!(j, Cga::e_plus());
    assert_eq!(eps, Cga::inf());
    assert_eq!(
        (scalar(2) + &j * 3_i64) * (scalar(4) + &j * 5_i64),
        scalar(23) + &j * 22_i64
    );
    assert_eq!(
        (scalar(2) + &eps * 3_i64) * (scalar(4) + &eps * 5_i64),
        scalar(8) + &eps * 22_i64
    );
    assert_eq!((scalar(1) + &j) * (scalar(1) - &j), scalar(0));
    assert_eq!(&eps * &eps, scalar(0));
    assert!(!eps.is_zero());
    assert_eq!(j.clifford_conjugate(), -&j);
    assert_eq!(eps.clifford_conjugate(), -&eps);
}

#[test]
fn mixed_unit_commutation_relations() {
    let i = Cga::ci();
    let i3 = Cga::pseudo3d();
    let j = Cga::sj();
    let eps = Cga::eps();
    let origin = Cga::orig();
    for element in [
        Cga::e1(),
        Cga::e2(),
        Cga::e3(),
        Cga::e_plus(),
        Cga::e_minus(),
        Cga::qi(),
        Cga::qj(),
        Cga::qk(),
        i3.clone(),
        j.clone(),
        eps.clone(),
        origin.clone(),
    ] {
        assert_eq!(&i * &element, &element * &i);
    }
    for vector in [Cga::e1(), Cga::e2(), Cga::e3()] {
        assert_eq!(&i3 * &vector, &vector * &i3);
        assert_eq!(&j * &vector, -(&vector * &j));
        assert_eq!(&eps * &vector, -(&vector * &eps));
        assert_eq!(&origin * &vector, -(&vector * &origin));
    }
    for quaternion in [Cga::qi(), Cga::qj(), Cga::qk()] {
        assert_eq!(&j * &quaternion, &quaternion * &j);
        assert_eq!(&eps * &quaternion, &quaternion * &eps);
        assert_eq!(&origin * &quaternion, &quaternion * &origin);
    }
    assert_eq!(&i3 * &j, -(&j * &i3));
    assert_eq!(&i3 * &eps, -(&eps * &i3));
    assert_eq!(&i3 * &origin, -(&origin * &i3));
    assert_eq!(&j * &eps + &eps * &j, scalar(2));
    assert_eq!(&j * &origin + &origin * &j, scalar(-1));
    assert_eq!(&origin * &eps + &eps * &origin, scalar(-2));
    assert_ne!(&j * &eps, &eps * &j);
}

#[test]
fn reversion_and_conjugation_have_distinct_meanings() {
    let i = Cga::ci();
    let i3 = Cga::pseudo3d();
    assert_eq!(i.reverse(), i);
    assert_eq!(i.grade_involution(), -&i);
    assert_eq!(i.clifford_conjugate(), -&i);
    assert_eq!(i3.reverse(), -&i3);
    assert_eq!(i3.clifford_conjugate(), i3);
    let z = scalar(2) + &i * 3_i64;
    assert_eq!(z.scalar_product(&z.reverse()), Real::from_int(-5));
    assert_eq!(z.norm_sq(), Real::from_int(13));
    assert_eq!(&z * &z.clifford_conjugate(), scalar(13));
}

#[test]
fn conjugate_norm_matches_each_catalog_number_system() {
    for (value, expected) in [
        (scalar(2) + Cga::ci() * 3_i64, 13),
        (
            scalar(2) + Cga::qi() * 3_i64 + Cga::qj() * 4_i64 + Cga::qk() * 5_i64,
            54,
        ),
        (scalar(2) + Cga::sj() * 3_i64, -5),
        (scalar(2) + Cga::eps() * 3_i64, 4),
        (Cga::ci(), 1),
        (Cga::qi(), 1),
        (Cga::sj(), -1),
        (Cga::eps(), 0),
    ] {
        assert_eq!(value.norm_sq(), Real::from_int(expected));
        let wrapped: Number = value.into();
        assert_eq!(
            wrapped.norm_sq().extract_complex().expect("scalar norm"),
            Number::from_int(expected)
                .extract_complex()
                .expect("real norm")
        );
    }
    assert_eq!((scalar(1) + Cga::sj()).norm_sq(), Real::zero());
}

#[test]
fn complex_scalars_map_to_the_central_blade() {
    let coefficient = Number::from_complex(Number::from_int(2), Number::from_int(3));
    let value = CliffordNumber::scalar(Cga::gens(), coefficient.clone()).expect("central scalar");
    assert_eq!(value.norm_sq(), Real::from_int(13));
    assert_eq!(value.coeff(0), &Real::from_int(2));
    assert_eq!(value.coeff(31), &Real::from_int(3));
    assert_eq!(
        Number::from(value.clifford_conjugate())
            .extract_complex()
            .expect("complex scalar"),
        coefficient.conj().extract_complex().expect("complex value")
    );
    let mixed = Number::from(Cga::e_plus()) * &coefficient;
    let expected = Number::from(
        Cga::e_plus() * Real::from_int(2) + Cga::e_plus() * Cga::ci() * Real::from_int(3),
    );
    assert_eq!(mixed, expected);
    assert_eq!(Number::from(Cga::ci()), Number::i());
}

#[test]
fn multivector_projections_do_not_discard_other_blades_during_extraction() {
    let value: Number = (scalar(2) + Cga::ci() * 3_i64 + Cga::e1() * 4_i64).into();
    assert_eq!(value.re(), Number::from_int(2));
    assert_eq!(value.im(), Number::from_int(3));
    assert!(value.is_int_ring());
    assert!(!value.is_int());
    assert!(!value.is_real());
    assert!(!value.is_complex());
    assert_eq!(
        value.extract_complex(),
        Err(crate::NumAnafisError::ExpectedScalar)
    );
    assert!(
        !Number::from(Cga::e1() * Real::from_rational(crate::RationalType::new(3, 2)))
            .is_int_ring()
    );
    for unit in [Cga::qi(), Cga::sj(), Cga::eps()] {
        let wrapped: Number = unit.into();
        assert!(wrapped.im().is_zero());
        assert!(wrapped.is_int_ring());
    }
    let unit: Number = Cga::ci().into();
    assert_eq!(unit.im(), Number::one());
    assert!(unit.is_complex());
    assert!(unit.is_int_ring());
    assert_eq!(
        unit.extract_complex().expect("complex unit").im,
        Real::one()
    );
    let gens =
        crate::GeneratorSet::new(vec![(0, 1)]).expect("valid generator labels, metrics, and count");
    let no_complex_unit: Number = CliffordNumber::generator(gens, 0)
        .expect("generator")
        .into();
    assert!(no_complex_unit.im().is_zero());
    assert!(!no_complex_unit.is_real());
}

#[test]
fn clifford_extrema_compare_only_real_scalars() {
    assert_eq!(scalar(-3).max(&scalar(-2)), scalar(-2));
    assert_eq!(scalar(3).min(&scalar(2)), scalar(2));
    for value in [Cga::ci(), Cga::e1()] {
        assert!(!value.max(&scalar(0)).is_finite());
        assert!(!scalar(0).min(&value).is_finite());
        let wrapped: Number = value.into();
        assert!(wrapped.max(&Number::zero()).is_nan());
        assert!(Number::zero().min(&wrapped).is_nan());
    }
}

#[test]
fn all_cga_blade_products_match_generator_word_reduction() {
    let blades: Vec<_> = (0..32).map(blade).collect();
    for (a, lhs) in blades.iter().enumerate() {
        for (b, rhs) in blades.iter().enumerate() {
            let (mask, sign) = reduce_word(a, b);
            let expected = &blades[mask] * sign;
            assert_eq!(lhs * rhs, expected, "blade product {a} * {b}");
        }
    }
}

proptest! {
    #[test]
    fn conformal_points_are_null_and_encode_distance(
        x in prop::array::uniform3(-10_i64..=10),
        y in prop::array::uniform3(-10_i64..=10),
    ) {
        let px = point(x);
        let py = point(y);
        let distance_squared: i64 = x.iter().zip(y).map(|(a, b)| (a-b)*(a-b)).sum();
        prop_assert_eq!(&px * &px, scalar(0));
        prop_assert_eq!(px.scalar_product(&Cga::inf()), Real::from_int(-1));
        prop_assert_eq!(&px * &py + &py * &px, scalar(-distance_squared));
    }

    #[test]
    fn catalog_linear_combinations_obey_associativity_and_distributivity(
        a in prop::array::uniform3(-3_i64..=3),
        b in prop::array::uniform3(-3_i64..=3),
        c in prop::array::uniform3(-3_i64..=3),
    ) {
        let lhs = scalar(a[0]) + Cga::ci() * a[1] + Cga::qi() * a[2];
        let middle = scalar(b[0]) + Cga::sj() * b[1] + Cga::qj() * b[2];
        let rhs = scalar(c[0]) + Cga::eps() * c[1] + Cga::qk() * c[2];
        prop_assert_eq!((&lhs * &middle) * &rhs, &lhs * (&middle * &rhs));
        prop_assert_eq!(&lhs * (&middle + &rhs), &lhs * &middle + &lhs * &rhs);
        prop_assert_eq!((&lhs * &middle).reverse(), middle.reverse() * lhs.reverse());
        prop_assert_eq!((&lhs * &middle).clifford_conjugate(), middle.clifford_conjugate() * lhs.clifford_conjugate());
    }

    #[test]
    fn subalgebra_norms_match_exact_formulas(a in -20_i64..=20, b in -20_i64..=20,
        c in -20_i64..=20, d in -20_i64..=20) {
        let complex = scalar(a) + Cga::ci() * b;
        let split = scalar(a) + Cga::sj() * b;
        let dual = scalar(a) + Cga::eps() * b;
        let quaternion = scalar(a) + Cga::qi() * b + Cga::qj() * c + Cga::qk() * d;
        prop_assert_eq!(complex.norm_sq(), Real::from_int(a*a+b*b));
        prop_assert_eq!(split.norm_sq(), Real::from_int(a*a-b*b));
        prop_assert_eq!(dual.norm_sq(), Real::from_int(a*a));
        prop_assert_eq!(quaternion.norm_sq(), Real::from_int(a*a+b*b+c*c+d*d));
    }
}

fn point(x: [i64; 3]) -> CliffordNumber {
    let norm: i64 = x.iter().map(|value| value * value).sum();
    let vector = Cga::e1() * x[0] + Cga::e2() * x[1] + Cga::e3() * x[2];
    vector
        + Cga::orig()
        + Cga::inf() * Real::from_rational(crate::types::RationalType::new(norm, 2))
}

fn scalar(value: i64) -> CliffordNumber {
    CliffordNumber::scalar(Cga::gens(), value).expect("valid scalar")
}

fn blade(mask: usize) -> CliffordNumber {
    let mut value = scalar(0);
    value.set_coeff(mask, 1);
    value
}

// Independent reference: sort a concatenated generator word by adjacent
// swaps, then contract identical adjacent generators using the metric.
fn reduce_word(a: usize, b: usize) -> (usize, i64) {
    let mut word = vec![];
    for mask in [a, b] {
        word.extend((0..5).filter(|position| mask & (1 << position) != 0));
    }
    let mut sign = 1_i64;
    for end in (1..word.len()).rev() {
        for index in 0..end {
            if word[index] > word[index + 1] {
                word.swap(index, index + 1);
                sign = -sign;
            }
        }
    }
    let mut mask = 0;
    let mut index = 0;
    while index < word.len() {
        if word.get(index + 1) == Some(&word[index]) {
            if word[index] == 4 {
                sign = -sign;
            }
            index += 2;
        } else {
            mask |= 1 << word[index];
            index += 1;
        }
    }
    (mask, sign)
}

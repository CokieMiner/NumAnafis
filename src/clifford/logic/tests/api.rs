//! Constructor validation, numerical ordering, and contraction contracts.

use core::cmp::Ordering;

use alloc::{vec, vec::Vec};

use proptest::prelude::*;

use crate::{Cga, CliffordNumber, FastClifford, GeneratorSet, NumAnafisError, Number, Real};

#[test]
fn generator_construction_sorts_and_validates() {
    let gens = GeneratorSet::new(vec![(302, -1), (300, 0), (301, 1)])
        .expect("distinct labels and valid metrics");
    assert_eq!((gens.id_at(0), gens.metric_at(0)), (300, 0));
    assert_eq!((gens.id_at(1), gens.metric_at(1)), (301, 1));
    assert_eq!((gens.id_at(2), gens.metric_at(2)), (302, -1));
    assert!(matches!(
        GeneratorSet::new(vec![(300, 1), (300, -1)]),
        Err(NumAnafisError::DuplicateGenerator { id: 300 })
    ));
    assert!(matches!(
        GeneratorSet::new(vec![(300, 2)]),
        Err(NumAnafisError::InvalidGeneratorMetric { id: 300, metric: 2 })
    ));
    assert!(matches!(
        GeneratorSet::new((0..usize::BITS).map(|id| (id, 1)).collect()),
        Err(NumAnafisError::ActiveGeneratorsTooLargeForPlatform { .. })
    ));
}

#[test]
#[should_panic(expected = "union blade count exceeds platform capacity")]
fn generator_union_preserves_platform_capacity() {
    let midpoint = usize::BITS >> 1;
    let left = GeneratorSet::new((0..midpoint).map(|id| (id, 1)).collect())
        .expect("half the index width is representable");
    let right = GeneratorSet::new((midpoint..usize::BITS).map(|id| (id, 1)).collect())
        .expect("half the index width is representable");
    let _union = left.union_with(&right);
}

#[test]
fn dense_coefficients_use_active_counts_and_preserve_representations() {
    let gens = GeneratorSet::new(vec![(300, 1)]).expect("valid algebra");
    let value =
        CliffordNumber::from_coeffs(gens.clone(), [Real::from_int(3), Real::from_float(-0.0)])
            .expect("two active coefficients");
    assert_eq!(value.coeffs_slice().len(), 2);
    assert_eq!(*value.coeff(0), Real::from_int(3));
    assert!(matches!(value.coeff(1), Real::Float(component) if component.is_sign_negative()));
    for coefficients in [vec![], vec![Real::zero()], vec![Real::zero(); 3]] {
        let found = coefficients.len();
        let result = CliffordNumber::from_coeffs(gens.clone(), coefficients);
        assert!(
            matches!(result, Err(NumAnafisError::DenseCoefficientLengthMismatch(details))
            if details.expected == 2 && details.found == found)
        );
    }
    let mut intermittent = [Some(Real::one()), None, Some(Real::one())].into_iter();
    let intermittent_coefficients = core::iter::from_fn(move || intermittent.next().flatten());
    assert!(matches!(
        CliffordNumber::from_coeffs(gens, intermittent_coefficients),
        Err(NumAnafisError::DenseCoefficientLengthMismatch(details))
            if details.expected == 2 && details.found == 1
    ));
    let heap_gens = GeneratorSet::new((300..306).map(|id| (id, 1)).collect())
        .expect("six generators fit the index width");
    let coefficients: Vec<Real> = (0..64).map(Real::from_int).collect();
    let heap_value = CliffordNumber::from_coeffs(heap_gens, coefficients.clone())
        .expect("64 active coefficients");
    assert_eq!(heap_value.coeffs_slice(), coefficients);
    let scalar = CliffordNumber::from_coeffs(GeneratorSet::empty(), [Real::one()])
        .expect("empty algebra has one blade");
    assert!(scalar.is_one());
}

#[test]
fn sparse_terms_accumulate_and_validate_masks() {
    let gens = GeneratorSet::new(vec![(300, 1)]).expect("valid algebra");
    let value = CliffordNumber::from_sparse(
        gens.clone(),
        [
            (1, Real::from_int(2)),
            (0, Real::from_int(3)),
            (1, Real::from_int(-2)),
        ],
    )
    .expect("valid masks");
    assert_eq!(*value.coeff(0), Real::from_int(3));
    assert!(value.coeff(1).is_zero());
    assert!(
        CliffordNumber::from_sparse(gens.clone(), [])
            .expect("empty sum")
            .is_zero()
    );
    assert!(matches!(
        CliffordNumber::from_sparse(gens, [(2, Real::one())]),
        Err(NumAnafisError::BladeIndexOutOfRange { blade: 2, count: 2 })
    ));
}

#[test]
fn numerical_order_and_extrema_reject_non_real_values() {
    let gens = Cga::gens();
    let small = CliffordNumber::scalar(gens.clone(), 2).expect("real scalar");
    let large = CliffordNumber::scalar(gens.clone(), 3).expect("real scalar");
    assert_eq!(small.partial_cmp(&large), Some(Ordering::Less));
    assert_eq!(small.max(&large), large);
    assert_eq!(small.min(&large), small);
    for unordered in [Cga::e1(), Cga::ci(), CliffordNumber::nan(gens.clone())] {
        let reflexive_order = (!unordered.coeff(0).is_nan()).then_some(Ordering::Equal);
        assert_eq!(unordered.partial_cmp(&unordered), reflexive_order);
        assert_eq!(unordered.partial_cmp(&small), None);
        assert_eq!(small.partial_cmp(&unordered), None);
        for result in [
            unordered.max(&small),
            unordered.min(&small),
            small.max(&unordered),
            small.min(&unordered),
        ] {
            assert!(result.coeff(0).is_nan());
            assert_eq!(result.generator_set(), &gens);
        }
        let number = Number::from(unordered);
        assert_eq!(number.partial_cmp(&number), reflexive_order);
        assert!(number.max(&number).is_nan());
        assert!(number.min(&number).is_nan());
        assert_eq!(number.partial_cmp(&Number::from(2)), None);
        assert!(number.max(&Number::from(2)).is_nan());
        assert!(number.min(&Number::from(2)).is_nan());
    }
    let other_context = CliffordNumber::scalar(GeneratorSet::empty(), 2).expect("real scalar");
    assert_eq!(small.partial_cmp(&other_context), None);
    assert_eq!(
        Number::from(small).partial_cmp(&Number::from(other_context)),
        Some(Ordering::Equal)
    );
    assert_eq!(Cga::e1().total_cmp(&Cga::e2()), Ordering::Greater);
}

#[test]
fn fixed_signature_order_is_real_only() {
    let mut small = FastClifford::<3, 0, 0>::zero();
    let mut large = small.clone();
    small.coeffs_mut_slice()[0] = Real::from_int(2);
    large.coeffs_mut_slice()[0] = Real::from_int(3);
    assert_eq!(small.partial_cmp(&large), Some(Ordering::Less));
    large.coeffs_mut_slice()[1] = Real::one();
    assert_eq!(large.partial_cmp(&large), Some(Ordering::Equal));
    assert_eq!(small.partial_cmp(&large), None);
    assert_eq!(small.total_cmp(&large), Ordering::Less);
}

#[test]
fn contraction_direction_and_scalar_part_are_distinct() {
    let vector = Cga::e1();
    let bivector = &vector * &Cga::e2();
    assert_eq!(vector.left_contraction(&bivector), Cga::e2());
    assert!(bivector.left_contraction(&vector).is_zero());
    assert_eq!(vector.scalar_product(&bivector), Real::zero());
    assert_eq!(Cga::e_minus().left_contraction(&Cga::e_minus()), scalar(-1));
    assert_eq!(scalar(2).left_contraction(&vector), &vector * 2_i64);
    assert!(vector.left_contraction(&scalar(2)).is_zero());
}

proptest! {
    #[test]
    fn sparse_and_dense_construction_agree(coefficients in prop::collection::vec(-100_i64..100, 8)) {
        let gens = GeneratorSet::new(vec![(300, 1), (301, -1), (302, 0)]).expect("valid algebra");
        let sparse = CliffordNumber::from_sparse(gens.clone(), coefficients.iter().enumerate().map(|(mask, value)| (mask, Real::from_int(*value)))).expect("valid masks");
        let dense = CliffordNumber::from_coeffs(gens, coefficients.into_iter().map(Real::from_int)).expect("eight coefficients");
        prop_assert_eq!(sparse, dense);
    }

    #[test]
    fn contraction_matches_homogeneous_grade_projection(left_mask in 0_usize..8, right_mask in 0_usize..8, left_coefficient in -10_i64..10, right_coefficient in -10_i64..10) {
        let gens = GeneratorSet::new(vec![(300, 1), (301, -1), (302, 0)]).expect("valid algebra");
        let left = CliffordNumber::from_sparse(gens.clone(), [(left_mask, Real::from_int(left_coefficient))]).expect("valid mask");
        let right = CliffordNumber::from_sparse(gens.clone(), [(right_mask, Real::from_int(right_coefficient))]).expect("valid mask");
        let expected = if left_mask.count_ones() <= right_mask.count_ones() {
            (&left * &right).grade(right_mask.count_ones() - left_mask.count_ones())
        } else {
            CliffordNumber::zero(gens).expect("valid algebra")
        };
        prop_assert_eq!(left.left_contraction(&right), expected);
    }
}

fn scalar(value: i64) -> CliffordNumber {
    CliffordNumber::scalar(Cga::gens(), value).expect("real scalar")
}

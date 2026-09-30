use core::{array, cmp::Ordering};

use alloc::{boxed::Box, vec, vec::Vec};

use crate::{
    error::{DenseLengthError, IndexOutOfRangeError, NumAnafisError},
    types::{Number, Real},
};

use super::{
    BladeAlgebra, CliffordCoeffs, CliffordNumber, GeneratorSet, INLINE_COEFF_COUNT,
    INLINE_GENERATOR_LIMIT,
};
// CliffordNumber — internal helpers

impl CliffordNumber {
    /// Scalar Cartesian parts in the declared central complex subalgebra.
    pub(crate) fn central_scalar_parts(&self) -> Option<(&Real, &Real)> {
        let imaginary = self.gens.complex_blade();
        if self
            .nonzero_blades()
            .any(|(blade, _)| blade != 0 && Some(blade) != imaginary)
        {
            return None;
        }
        Some((
            self.coeff(0),
            imaginary.map_or(&Real::Int(0), |blade| self.coeff(blade)),
        ))
    }

    /// Borrows the active blade coefficients; padding is excluded.
    #[must_use]
    pub fn coeffs_slice(&self) -> &[Real] {
        let count = self.blade_count();
        match self.coeffs {
            CliffordCoeffs::Inline(ref a) => &a[..count],
            CliffordCoeffs::Heap(ref v) => v,
        }
    }

    /// Mutably borrows the active blade coefficients; padding is excluded.
    #[must_use]
    pub fn coeffs_mut_slice(&mut self) -> &mut [Real] {
        let count = self.blade_count();
        match self.coeffs {
            CliffordCoeffs::Inline(ref mut a) => &mut a[..count],
            CliffordCoeffs::Heap(ref mut v) => v.as_mut_slice(),
        }
    }

    /// Initializes inline coefficients or allocates dense heap storage for n generators.
    #[must_use]
    fn make_zero_coeffs(n: u8) -> CliffordCoeffs {
        let z = Real::zero();
        if n <= INLINE_GENERATOR_LIMIT {
            CliffordCoeffs::Inline(array::from_fn(|_| z.clone()))
        } else {
            CliffordCoeffs::Heap(vec![z; BladeAlgebra::coeff_count_unchecked(n)])
        }
    }

    /// Zero multivector for `gens` without validating the generator count.
    ///
    /// Callers must ensure the blade count fits the coefficient storage.
    ///
    /// # Panics
    /// Panics if the generator count exceeds 255 or if `2^n` overflows `usize`.
    #[must_use]
    pub(crate) fn zero_unchecked(gens: GeneratorSet) -> Self {
        let n = u8::try_from(gens.len()).expect("active generator count <= 64 fits in u8");
        Self {
            coeffs: Self::make_zero_coeffs(n),
            gens,
        }
    }

    /// Re-embed `self` into a new generator set using a precomputed bit-permutation.
    ///
    /// # Panics
    /// Panics if either generator count is unrepresentable, or if `perm` maps
    /// a blade outside the new coefficient storage.
    #[must_use]
    pub(crate) fn reembed(&self, new_gens: GeneratorSet, perm: &[u8]) -> Self {
        let new_n = u8::try_from(new_gens.len()).expect("active generator count <= 64 fits in u8");
        let new_count = BladeAlgebra::coeff_count_unchecked(new_n);
        let old_count = self.blade_count();
        let z = Real::zero();

        let mut buf: Vec<Real> = (0..new_count).map(|_| z.clone()).collect();
        for old_blade in 0..old_count {
            let c = &self.coeffs_slice()[old_blade];
            if c.is_zero() {
                continue;
            }
            buf[BladeAlgebra::permute_blade(old_blade, perm)] = c.clone();
        }

        if new_n <= INLINE_GENERATOR_LIMIT {
            let arr: [Real; INLINE_COEFF_COUNT] = array::from_fn(|i| {
                if i < new_count {
                    buf[i].clone()
                } else {
                    z.clone()
                }
            });
            Self {
                coeffs: CliffordCoeffs::Inline(arr),
                gens: new_gens,
            }
        } else {
            Self {
                coeffs: CliffordCoeffs::Heap(buf),
                gens: new_gens,
            }
        }
    }
}
// CliffordNumber — public API

impl CliffordNumber {
    // --- Inspection ---

    /// Number of basis blades = 2^n for the active generator set.
    ///
    /// This is the count of *valid* blades, not the backing storage length:
    /// inline storage always holds 32 coefficients regardless of `n`.
    #[must_use]
    pub const fn blade_count(&self) -> usize {
        1_usize << self.gens.len()
    }

    /// Number of active generators.
    #[must_use]
    pub const fn n_generators(&self) -> usize {
        self.gens.len()
    }

    /// The active [`GeneratorSet`].
    #[must_use]
    pub const fn generator_set(&self) -> &GeneratorSet {
        &self.gens
    }

    /// Coefficient at blade bitmask `blade`.
    ///
    /// # Panics
    /// Panics if `blade` is out of range for this algebra (`blade >= 2^n`,
    /// including inline padding lanes).
    #[must_use]
    pub fn coeff(&self, blade: usize) -> &Real {
        assert!(
            blade < self.blade_count(),
            "blade index {blade} out of range for {}-generator algebra",
            self.gens.len()
        );
        &self.coeffs_slice()[blade]
    }

    /// Set the coefficient at blade bitmask `blade`.
    ///
    /// # Panics
    /// Panics if `blade` is out of range for this algebra (`blade >= 2^n`,
    /// including inline padding lanes).
    ///
    /// ```compile_fail
    /// use num_anafis::{Cga, i};
    /// let mut value = Cga::e1();
    /// value.set_coeff(1, i()); // A complex scalar is not a real coefficient.
    /// ```
    pub fn set_coeff(&mut self, blade: usize, value: impl Into<Real>) {
        assert!(
            blade < self.blade_count(),
            "blade index {blade} out of range for {}-generator algebra",
            self.gens.len()
        );
        self.coeffs_mut_slice()[blade] = value.into();
    }

    /// Iterator over `(blade_mask, &coefficient)` for every non-zero blade.
    pub fn nonzero_blades(&self) -> impl Iterator<Item = (usize, &Real)> {
        self.coeffs_slice()
            .iter()
            .take(self.blade_count())
            .enumerate()
            .filter(|&(_, c)| !c.is_zero())
    }

    // --- Constructors ---

    /// All-zero multivector in `gens`.
    ///
    /// # Errors
    /// [`NumAnafisError::ActiveGeneratorsTooLargeForPlatform`] if 2ⁿ overflows `usize`
    /// or the generator count exceeds `u8::MAX`.
    pub fn zero(gens: GeneratorSet) -> Result<Self, NumAnafisError> {
        let n = gens.len();
        let n_u8 = u8::try_from(n).map_err(|_e| {
            NumAnafisError::ActiveGeneratorsTooLargeForPlatform { active: u8::MAX }
        })?;
        let _ = BladeAlgebra::coeff_count(n_u8)?;
        Ok(Self::zero_unchecked(gens))
    }

    /// Embeds an ordinary scalar, identifying its imaginary unit with the
    /// declared central pseudoscalar. Real inputs occupy grade zero.
    ///
    /// # Errors
    /// Same as [`Self::zero`], plus [`NumAnafisError::ExpectedScalar`] for a
    /// general multivector and [`NumAnafisError::UnsupportedComplexStructure`]
    /// for a complex input in an algebra without a compatible central unit.
    pub fn scalar(gens: GeneratorSet, value: impl Into<Number>) -> Result<Self, NumAnafisError> {
        let mut mv = Self::zero(gens)?;
        let parts = value.into().into_complex()?;
        if !parts.im.is_zero()
            || matches!(&parts.im,Real::Float(imaginary) if imaginary.is_sign_negative())
        {
            let blade = mv
                .gens
                .complex_blade()
                .ok_or(NumAnafisError::UnsupportedComplexStructure)?;
            mv.coeffs_mut_slice()[blade] = parts.im;
        }
        mv.coeffs_mut_slice()[0] = parts.re;
        Ok(mv)
    }

    /// NaN multivector in `gens` without validating the generator count.
    ///
    /// Scalar blade is NaN, all other blades are zero. Callers must ensure
    /// the blade count fits the coefficient storage.
    ///
    /// # Panics
    /// Panics if the generator count exceeds 255 or if `2^n` overflows `usize`.
    #[must_use]
    pub fn nan(gens: GeneratorSet) -> Self {
        let mut mv = Self::zero_unchecked(gens);
        mv.coeffs_mut_slice()[0] = Real::nan();
        mv
    }

    /// Basis generator vector with only bit `index` set.
    ///
    /// # Errors
    /// [`NumAnafisError::GeneratorIndexOutOfRange`] when `index >= n`.
    /// [`NumAnafisError::ActiveGeneratorsTooLargeForPlatform`] if the generator count
    /// overflows `u8` or 2ⁿ overflows `usize`.
    pub fn generator(gens: GeneratorSet, index: u8) -> Result<Self, NumAnafisError> {
        let n = gens.len();
        let n_u8 = u8::try_from(n).map_err(|_e| {
            NumAnafisError::ActiveGeneratorsTooLargeForPlatform { active: u8::MAX }
        })?;
        if usize::from(index) >= n {
            return Err(NumAnafisError::GeneratorIndexOutOfRange(Box::new(
                IndexOutOfRangeError {
                    active: n_u8,
                    index,
                },
            )));
        }
        let mut mv = Self::zero(gens)?;
        mv.coeffs_mut_slice()[1_usize << usize::from(index)] = Real::one();
        Ok(mv)
    }

    /// Builds from exactly `2^n` real coefficients in blade-mask order.
    /// Arrays and vectors transfer their coefficients without cloning.
    ///
    /// # Errors
    /// [`NumAnafisError::DenseCoefficientLengthMismatch`] for an incorrect count.
    pub fn from_coeffs(
        gens: GeneratorSet,
        coeffs: impl IntoIterator<Item = Real>,
    ) -> Result<Self, NumAnafisError> {
        let expected = 1_usize << gens.len();
        let mut coefficients = coeffs.into_iter().fuse();
        let (storage, found) = if gens.len() <= usize::from(INLINE_GENERATOR_LIMIT) {
            let mut inline = array::from_fn(|_| Real::zero());
            let mut found = 0;
            for (target, value) in inline.iter_mut().take(expected).zip(coefficients.by_ref()) {
                *target = value;
                found += 1;
            }
            (CliffordCoeffs::Inline(inline), found + coefficients.count())
        } else {
            let dense: Vec<Real> = coefficients.collect();
            let found = dense.len();
            (CliffordCoeffs::Heap(dense), found)
        };
        if found != expected {
            return Err(NumAnafisError::DenseCoefficientLengthMismatch(Box::new(
                DenseLengthError { expected, found },
            )));
        }
        Ok(Self {
            coeffs: storage,
            gens,
        })
    }

    /// Builds from `(blade_mask, coefficient)` terms; repeated masks are added.
    /// Unspecified blades have zero coefficients.
    ///
    /// # Errors
    /// [`NumAnafisError::BladeIndexOutOfRange`] for a mask outside `0..2^n`.
    pub fn from_sparse(
        gens: GeneratorSet,
        terms: impl IntoIterator<Item = (usize, Real)>,
    ) -> Result<Self, NumAnafisError> {
        let mut result = Self::zero_unchecked(gens);
        let count = result.blade_count();
        for (blade, coefficient) in terms {
            if blade >= count {
                return Err(NumAnafisError::BladeIndexOutOfRange { blade, count });
            }
            result.coeffs_mut_slice()[blade] += &coefficient;
        }
        Ok(result)
    }
}

// Strict numerical ordering requires real scalars in the same algebra.

impl PartialOrd for CliffordNumber {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.gens != other.gens {
            return None;
        }
        match (self.central_scalar_parts(), other.central_scalar_parts()) {
            (Some((a, b)), Some((c, d))) if b.is_zero() && d.is_zero() => a.partial_cmp(c),
            _ => (self == other).then_some(Ordering::Equal),
        }
    }
}

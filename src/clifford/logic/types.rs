use core::{
    cmp::Ordering,
    hash::{Hash, Hasher},
};

use alloc::{
    string::{String, ToString},
    sync::Arc,
    vec,
    vec::Vec,
};

use crate::{error::NumAnafisError, types::Real};

use super::{E_MINUS, E_PLUS, E1, E2, E3};
// Storage and precomputation limits

/// Maximum generators for inline (stack-resident) coefficient storage.
pub const INLINE_GENERATOR_LIMIT: u8 = 5;
/// Inline coefficient count = 2^5 = 32.
pub const INLINE_COEFF_COUNT: usize = 32;
/// Maximum generator count with a precomputed Cayley sign table.
///
/// The table has at most 128^2 entries; larger signatures compute signs directly.
const CAYLEY_GENERATOR_LIMIT: usize = 7;
// GeneratorSet

/// A sorted, duplicate-free set of `(generator_id, metric)` pairs defining an algebra.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct GeneratorSet {
    entries: Vec<(u32, i8)>,
    /// Precomputed Cayley sign table for blade products.
    ///
    /// Empty when the generator count exceeds [`CAYLEY_GENERATOR_LIMIT`];
    /// geometric products then fall back to [`mul_blades`] directly.
    pub(crate) cayley_signs: Arc<[i8]>,
}

impl PartialEq for GeneratorSet {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries
    }
}

impl Eq for GeneratorSet {}

impl Hash for GeneratorSet {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.entries.hash(state);
    }
}

fn compute_cayley_signs(entries: &[(u32, i8)]) -> Arc<[i8]> {
    let n = entries.len();
    if n > CAYLEY_GENERATOR_LIMIT {
        return Arc::from(Vec::new());
    }
    let count = 1 << n;
    let mut table = vec![0; count * count];
    let dummy = GeneratorSet {
        entries: entries.to_vec(),
        cayley_signs: Arc::from(Vec::new()),
    };
    for a in 0..count {
        for b in 0..count {
            if let Some((_, factor)) = BladeAlgebra::mul_blades(&dummy, a, b) {
                table[a * count + b] = factor;
            }
        }
    }
    Arc::from(table)
}

impl GeneratorSet {
    /// Oriented full pseudoscalar when it is central and squares to minus one.
    /// No degenerate or even-dimensional algebra has this declared scalar unit.
    #[must_use]
    pub fn complex_blade(&self) -> Option<usize> {
        let n = self.len();
        if n.is_multiple_of(2) || u32::try_from(n).map_or(true, |count| count >= usize::BITS) {
            return None;
        }
        let blade = (1_usize << n) - 1;
        match BladeAlgebra::mul_blades(self, blade, blade) {
            Some((0, -1)) => Some(blade),
            _ => None,
        }
    }

    /// Empty set — represents the pure-scalar subalgebra.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            entries: Vec::new(),
            cayley_signs: compute_cayley_signs(&[]),
        }
    }

    /// Builds an algebra from `(id, metric)` pairs, sorting by generator label.
    ///
    /// # Errors
    /// Rejects duplicate labels, metrics outside `{-1, 0, 1}`, and counts
    /// for which `2^n` cannot be represented by `usize`.
    pub fn new(mut entries: Vec<(u32, i8)>) -> Result<Self, NumAnafisError> {
        if u32::try_from(entries.len()).map_or(true, |count| count >= usize::BITS) {
            return Err(NumAnafisError::ActiveGeneratorsTooLargeForPlatform {
                active: u8::try_from(entries.len()).unwrap_or(u8::MAX),
            });
        }
        for &(id, metric) in &entries {
            if !matches!(metric, -1..=1) {
                return Err(NumAnafisError::InvalidGeneratorMetric { id, metric });
            }
        }
        entries.sort_unstable_by_key(|&(id, _)| id);
        if let Some(pair) = entries.windows(2).find(|pair| pair[0].0 == pair[1].0) {
            return Err(NumAnafisError::DuplicateGenerator { id: pair[0].0 });
        }
        Ok(Self::from_validated(entries))
    }

    /// Constructs from strictly increasing labels, valid metrics, and a count
    /// below `usize::BITS`; callers establish these invariants structurally.
    pub(crate) fn from_validated(entries: Vec<(u32, i8)>) -> Self {
        Self {
            cayley_signs: compute_cayley_signs(&entries),
            entries,
        }
    }

    /// Number of active generators.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.entries.len()
    }

    /// `true` when there are no generators (scalar subalgebra).
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Generator ID at bit position `i`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    #[must_use]
    pub fn id_at(&self, i: usize) -> u32 {
        self.entries[i].0
    }

    /// Metric factor (+1, −1, or 0) of the generator at bit position `i`.
    ///
    /// # Panics
    /// Panics if `i` is out of bounds.
    #[must_use]
    pub fn metric_at(&self, i: usize) -> i8 {
        self.entries[i].1
    }

    /// Compute the union of `self` and `other`.
    ///
    /// Returns `(union, perm_self, perm_other)` where `perm_x[i]` is the bit position of
    /// `x`'s i-th generator in the union.
    ///
    /// # Panics
    ///
    /// Panics if the same generator ID exists in both sets but with different metrics.
    /// Panics if `2^n` for the union cannot be represented by `usize`.
    #[must_use]
    pub fn union_with(&self, other: &Self) -> (Self, Vec<u8>, Vec<u8>) {
        let mut union = Vec::with_capacity(self.entries.len() + other.entries.len());
        let mut perm_a = vec![0_u8; self.entries.len()];
        let mut perm_b = vec![0_u8; other.entries.len()];

        let (mut ia, mut ib) = (0_usize, 0_usize);
        let mut pos: u8 = 0;

        while ia < self.entries.len() && ib < other.entries.len() {
            let (aid, am) = self.entries[ia];
            let (bid, bm) = other.entries[ib];

            match aid.cmp(&bid) {
                Ordering::Less => {
                    union.push((aid, am));
                    perm_a[ia] = pos;
                    ia += 1;
                }
                Ordering::Equal => {
                    assert_eq!(am, bm, "same generator ID must have the same metric");
                    union.push((aid, am));
                    perm_a[ia] = pos;
                    perm_b[ib] = pos;
                    ia += 1;
                    ib += 1;
                }
                Ordering::Greater => {
                    union.push((bid, bm));
                    perm_b[ib] = pos;
                    ib += 1;
                }
            }
            pos = pos.checked_add(1).expect("union exceeds 255 generators");
        }
        for (idx, &e) in self.entries[ia..].iter().enumerate() {
            union.push(e);
            perm_a[ia + idx] = pos;
            pos = pos.checked_add(1).expect("union exceeds 255 generators");
        }
        for (idx, &e) in other.entries[ib..].iter().enumerate() {
            union.push(e);
            perm_b[ib + idx] = pos;
            pos = pos.checked_add(1).expect("union exceeds 255 generators");
        }

        assert!(
            u32::try_from(union.len()).is_ok_and(|count| count < usize::BITS),
            "union blade count exceeds platform capacity"
        );
        (Self::from_validated(union), perm_a, perm_b)
    }
}
// CliffordCoeffs — dense coefficient storage

#[expect(
    clippy::large_enum_variant,
    reason = "Inline path is intentionally stack-resident for n ≤ 5."
)]
#[derive(Debug, Clone, PartialEq, Hash)]
pub enum CliffordCoeffs {
    Inline([Real; INLINE_COEFF_COUNT]),
    Heap(Vec<Real>),
}
// CliffordNumber

/// A dense multivector in a Clifford algebra with labeled generators.
#[derive(Debug, Clone, PartialEq, Hash)]
#[non_exhaustive]
pub struct CliffordNumber {
    /// Blade coefficient storage (inline or heap).
    pub(crate) coeffs: CliffordCoeffs,
    /// Active generator set defining the algebra.
    pub(crate) gens: GeneratorSet,
}
// BladeAlgebra — blade-index products, sizing, permutations, and labels

/// Blade-index algebra: Cayley products, coefficient sizing, generator
/// permutations, and human-readable blade labels.
pub struct BladeAlgebra;

impl BladeAlgebra {
    pub fn coeff_count(n: u8) -> Result<usize, NumAnafisError> {
        1_usize
            .checked_shl(u32::from(n))
            .ok_or(NumAnafisError::ActiveGeneratorsTooLargeForPlatform { active: n })
    }

    pub fn coeff_count_unchecked(n: u8) -> usize {
        Self::coeff_count(n).expect("generator count was already validated")
    }

    /// Maps a blade mask from source generator positions to union positions.
    pub fn permute_blade(blade_old: usize, perm: &[u8]) -> usize {
        let mut new = 0_usize;
        for (i, &p) in perm.iter().enumerate() {
            if (blade_old >> i) & 1 == 1 {
                new |= 1_usize << usize::from(p);
            }
        }
        new
    }

    /// Compute the result blade and sign factor for the geometric product of blades `a` × `b`.
    pub fn mul_blades(gens: &GeneratorSet, a: usize, b: usize) -> Option<(usize, i8)> {
        let n = gens.len();

        let mut swaps = 0_u32;
        for i in 0..n {
            if (b >> i) & 1 == 1 {
                swaps += (a >> (i + 1)).count_ones();
            }
        }
        let mut factor: i8 = if swaps & 1 == 0 { 1 } else { -1 };

        let overlap = a & b;
        if overlap != 0 {
            for i in 0..n {
                if (overlap >> i) & 1 == 1 {
                    match gens.metric_at(i) {
                        1 => {}
                        -1 => factor = -factor,
                        // Zero metric (nilpotent) yields no product, as does any
                        // value outside the validated +1/-1/0 invariant.
                        _ => return None,
                    }
                }
            }
        }
        Some((a ^ b, factor))
    }

    /// Human-readable name for a generator ID.
    pub const fn gen_name(id: u32) -> &'static str {
        match id {
            E1 => "e1",
            E2 => "e2",
            E3 => "e3",
            E_PLUS => "e+",
            E_MINUS => "e-",
            _ => "e?",
        }
    }

    /// Basis blade label for display.
    pub fn blade_label(gens: &GeneratorSet, blade: usize) -> String {
        let mut label = String::new();
        for i in 0..gens.len() {
            if (blade >> i) & 1 == 1 {
                let name = Self::gen_name(gens.id_at(i));
                if name == "e?" {
                    label.push_str(&gens.id_at(i).to_string());
                } else {
                    label.push_str(name);
                }
            }
        }
        label
    }
}

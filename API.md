# Public Rust API

This guide covers every item re-exported by `src/lib.rs`, including its public
methods and operator support. Import these items directly from `num_anafis`.
The implementation modules are private.

## Exports and features

| Export | Purpose | Feature |
|---|---|---|
| `Number` | One numeric value for CAS literals and evaluation | Always |
| `Numeric` | Elementary evaluator function trait | Always |
| `n`, `r`, `c`, `i` | Convenience constructors returning `Number` | Always |
| `Real` | Integer, rational, or approximate real component | Always |
| `Complex` | Typed complex value with two real components | Always |
| `IntType` | Native integer alias: `i64` | Always |
| `RationalType` | Reduced rational with `i64` numerator and `u64` denominator | Always |
| `FloatType` | Native approximate alias: `f64` | Always |
| `F32Ext`, `F64Ext` | Primitive floating-point functions and real special kernels | Always |
| `NumAnafisError` | Checked construction and conversion failures | Always |
| `GeneratorSet` | Labeled Clifford generator metrics | `clifford` |
| `CliffordNumber` | Dynamic real-coefficient multivector | `clifford` |
| `Cga` | Conformal geometric algebra unit catalog | `clifford` |
| `FastClifford<P,Q,R>` | Multivector with a fixed signature of at most five generators | `clifford` |

`std` is enabled by default. Without it, the library uses `core` and `alloc`; the
consumer supplies an allocator. `python` enables PyO3 bindings and `std`, but adds
no Rust re-exports to this list. Arithmetic uses the native aliases in the table.
`epsilon()` returns machine epsilon for the float representation.

## Using Number in a CAS

The numeric literal payload can be `Number` for every supported domain:

```rust
use num_anafis::{Number, Numeric, c, i, n, r};

#[derive(Clone)]
enum Expr {
    Literal(Number),
    Symbol(String),
}

let literal = Expr::Literal(r(2, 3));
if let Expr::Literal(value) = literal {
    assert_eq!(value.pow(&n(-3)), r(27, 8));
}
assert_eq!(r(-4, 9).sqrt(), c(0, r(2, 3)));
assert_eq!(5 * i(), c(0, 5));
```

CAS evaluation uses the following contracts:

- Import `Numeric` for elementary function calls. Methods borrow their input and
  return a new value of the same implementing type.
- `Number` is `Clone`, not `Copy`. Borrow stored operands (`&a + &b`) and clone
  when ownership is actually required.
- Use `n(2)` or `2.into()` for an exact integer. `n(2.0)` remains a float.
  Parsing every literal through `f64` would lose the exact representation first.
- `is_int()` checks whether the represented value is a real integer, including
  finite floats with no fractional part. It cannot prove that a preceding
  approximate calculation was exact. `is_int_ring()` checks all components.
  A symbolic expression remains an AST node; the crate does not decide when to
  evaluate that node.
- `pow(&exponent)` replaces `powf(exponent)`. Names also differ for `expm1()`
  versus `exp_m1()` and `log1p()` versus `ln_1p()`.
- Numeric comparisons are partial. Strict order applies to real values; equal
  non-real values compare `Equal`, and distinct non-real values are unordered.
  `total_cmp()` provides a sorting order, not a mathematical complex order.
- `Hash` exists, but `Eq` and `Ord` do not exist for `Number`. A CAS using
  hash-map/set keys or interning needs its own explicit key equality policy;
  `Number` alone does not satisfy those containers' `Eq` requirement.
- Primitive projection is explicit. Check `is_real()` before using `to_f64()`
  when discarding an imaginary component would be invalid.
- `Number` has no assignment operators such as `+=`. Use
  `accumulator = &accumulator + &term`.

Once special functions are implemented and integrated into `Numeric` and its
dispatch, their values can use this same literal payload. Supporting complex
kernels alone does not automatically integrate Clifford functional calculus;
that evaluator path must also support each new function. AST algebraic rules
remain the CAS's responsibility.

## Number construction and arithmetic

### Convenience constructors

| Call | Return and behavior |
|---|---|
| `n(value)` | `Number`; accepts any `Into<Number>` value |
| `r(numer: i64, denom: i64)` | `Number`; reduced rational, accepts either denominator sign; zero denominator returns NaN |
| `c(re, im)` | `Number`; both arguments implement `Into<Number>`; constructs `re + i*im`, flattening scalar complex arguments |
| `i()` | `Number`; exact ordinary imaginary unit |

Use `n(2.0)` for a float and `5 * i()` for an imaginary multiple. `c` panics for a general multivector argument
that cannot be interpreted as an ordinary scalar.

### Standard conversions

`Number::from(value)`, `value.into()`, and `n(value)` accept `i8`, `i16`, `i32`,
`i64`, `u8`, `u16`, `u32`, `f32`, `f64`, `Real`, and `Complex`. With `clifford`,
they also accept an owned `CliffordNumber`. Converting a `Number` to itself uses
the standard identity conversion.

`Number::try_from(u64_value)` and `Number::try_from(usize_value)` return
`Result<Number, core::num::TryFromIntError>` and reject values above `i64::MAX`.
They never silently approximate an oversized integer. These two types do not
implement `Into<Number>`; use the checked conversion before calling `n` if needed.
`RationalType` enters through `Number::from_rational` or `Real::from_rational`.

### Associated constructors

| Method | Return | Meaning |
|---|---|---|
| `from_int(IntType)` | `Number` | Exact integer |
| `from_rational(RationalType)` | `Number` | Exact rational; denominator one becomes integer |
| `from_float(FloatType)` | `Number` | Approximate value, preserving integral floats |
| `from_complex(Number, Number)` | `Number` | Same scalar construction rule as `c` |
| `zero()`, `one()`, `neg_one()` | `Number` | Exact identities |
| `nan()` | `Number` | IEEE NaN |
| `pi()`, `e()`, `epsilon()` | `Number` | Approximate native constants |
| `i()` | `Number` | Exact imaginary unit |

A complex result with imaginary positive zero normalizes to a real value.
Imaginary negative floating zero is preserved for branch-side information.
Approximate components are not relabeled as exact integers during normalization.

### Operators

`+`, `-`, `*`, `/`, and unary `-` operate on owned and borrowed `Number` values.
Binary arithmetic also supports native `i64` and `f64` operands on either side,
including borrowed `Number` operands. Other primitive types can be passed through
`n` or standard conversions first.

```rust
use num_anafis::{i, n, r};

let a = r(2, 3);
let b = n(4);
assert_eq!(&a + &b, r(14, 3));
assert_eq!(3 / &n(2), r(3, 2));
assert_eq!(0.5 * i(), num_anafis::c(0, 0.5));
```

Native integer arithmetic wraps in debug and release. Rational arithmetic stays
exact within capacity and panics when its reduced numerator/denominator cannot
fit. Integer division creates a rational when needed. Float arithmetic follows
the approximate tier; division by zero can produce infinity or NaN. There are no
implicit overflow promotions, remainder operators, or assignment operators on
`Number`.

## Number inspection and conversion

`is_int()` tests the represented value without a tolerance. Thus `n(2.0)` passes
but `n(2.0000000000000004)` does not. Rounding may already have produced an
integral float, so this predicate is not proof of an exact mathematical result.
The float stays in the float tier. When a CAS requires an exact integer literal,
check `is_int() && !re().is_float()`.

`is_int()` identifies ordinary integers. `is_int_ring()` together with
`is_complex()` identifies non-real Gaussian integers; together with
`!is_real() && !is_complex()` it identifies integer-coefficient multivectors
outside the scalar complex subalgebra. Integral floats participate by stored
value, with the same rounding limitation.

All methods below take `&self` unless shown as constructors above.

| Method | Return | Contract |
|---|---|---|
| `is_rational()`, `is_float()` | `bool` | Test the active ordinary real storage tier |
| `is_clifford()` | `bool` | Test multivector storage; always false without the feature |
| `is_real()` | `bool` | Ordinary real scalar meaning, including a compatible scalar-only Clifford value |
| `is_complex()` | `bool` | Ordinary scalar meaning with a nonzero imaginary component |
| `is_int()` | `bool` | Real integer value, including finite integral floats; does not establish computation exactness |
| `is_int_ring()` | `bool` | Every component is an integer; complex scalars are Gaussian integers, multivectors have integer blade coefficients |
| `is_zero()`, `is_one()`, `is_neg_one()` | `bool` | Numeric identity tests |
| `is_negative()`, `is_positive()` | `bool` | Strict sign of a real scalar; false for non-real/general multivectors |
| `is_nan()`, `is_infinite()`, `is_finite()` | `bool` | Inspect all components |
| `re()`, `im()` | `Number` | Scalar and declared complex-unit coefficients, including general multivectors; CGA uses I5. `im()` is zero without a declared unit. Other blades are excluded |
| `conj()` | `Number` | Ordinary complex conjugation or its compatible Clifford anti-involution |
| `norm_sq()` | `Number` | Scalar squared modulus; a signed conjugate form on general multivectors |
| `arg()` | `Number` | Principal scalar argument; NaN for a general multivector |
| `to_f64()`, `to_f32()` | `f64`, `f32` | Approximate real-component projection; panics for a general multivector |
| `to_int()` | `IntType` | Truncating real-component projection; panics for a general multivector |
| `to_float()` | `Number` | Convert every component to the approximate tier |
| `trunc()` | `Number` | Truncate scalar components toward zero; unsupported Clifford inputs give NaN |
| `extract_complex()` | `Result<Complex, NumAnafisError>` | Convert the whole value to typed complex components; rejects extra blades, discards Clifford context |
| `as_clifford()` | `Option<&CliffordNumber>` | Borrow multivector storage; method exists only with `clifford` |
| `total_cmp(&other)` | `Ordering` | Scalars first, by approximate magnitude then exact component tie-breaks; general multivectors by context/components |
| `max(&other)`, `min(&other)` | `Number` | Real scalar comparison; non-real or unordered operands return NaN |
| `clamp(&min, &max)` | `Number` | Applies lower then upper real bounds; non-real/unordered inputs return NaN; bound ordering is not validated |
| `approx_eq_number(&other, &tolerance)` | `bool` | Finite scalar modulus-distance comparison, or coefficient distance for matching Clifford operands/tolerance |

Real finite equality and hashing agree across tiers without rounding exact
values first. Thus `n(1) == n(1.0)`, while `r(1,10) != n(0.1)`. Signed zeros
compare equal. NaN is nonreflexive and unordered. `Number` implements `Clone`,
`Debug`, `Display`, `Hash`, `PartialEq`, and `PartialOrd`.

`to_int()` saturates out-of-range float inputs like a Rust float-to-`i64` cast;
NaN projects to zero. It is a projection, not checked integer conversion.
`to_f64()` and `to_int()` read only the real component even for a non-real scalar.

## Numeric: elementary evaluator functions

Implementations exist for `Number`, `Real`, `Complex`, and, with `clifford`,
`CliffordNumber` and `FastClifford`. Primitive `f32`/`f64` use their native methods
and extension traits rather than implementing `Numeric`.

Unary methods have signature `fn name(&self) -> Self`.

| Methods | Meaning |
|---|---|
| `sin`, `cos`, `tan` | Trigonometric functions in radians |
| `cot`, `sec`, `csc` | Reciprocal trigonometric functions |
| `asin`, `acos`, `atan` | Inverse trigonometric functions |
| `acot`, `asec`, `acsc` | `atan(1/x)`, `acos(1/x)`, `asin(1/x)` |
| `sinh`, `cosh`, `tanh` | Hyperbolic functions |
| `coth`, `sech`, `csch` | Reciprocal hyperbolic functions |
| `asinh`, `acosh`, `atanh` | Inverse hyperbolic functions |
| `acoth`, `asech`, `acsch` | Inverse hyperbolic functions applied to `1/x` |
| `exp`, `expm1`, `exp_neg` | `exp(x)`, stable `exp(x)-1`, `exp(-x)` |
| `ln`, `log1p` | Natural logarithm, stable `ln(1+x)` |
| `sqrt`, `cbrt` | Square root, cube root |
| `abs` | Real absolute value or scalar complex modulus |
| `signum` | Real sign, or complex direction `z/abs(z)` for nonzero `z` |
| `floor`, `ceil`, `round`, `fract` | Rounding/fractional operations; complex scalars act componentwise |
| `negate` | Arithmetic negation, also available through unary `-` |
| `sinc` | Unnormalized `sin(x)/x`, with `sinc(0)=1` |

Binary methods have signature `fn name(&self, other: &Self) -> Self`:

| Method | Meaning |
|---|---|
| `atan2(&x)` | `self` is `y`; quadrant-aware angle for real scalar arguments; non-real/general Clifford inputs return NaN |
| `log_base(&base)` | `ln(self)/ln(base)` |
| `pow(&exponent)` | Power with exact integer paths and principal logarithmic continuation when needed |

`Number` promotes real arguments into complex evaluation when the supported
function leaves the real domain. A typed `Real` result must remain real and can
therefore be NaN outside that domain. Use `Number` for automatic domain promotion.
Complex logarithms and general noninteger powers use principal conventions;
`cbrt` on the real axis retains the real cube root, including negative inputs.
`Complex` has an inherent `abs() -> Real`; `Numeric::abs(&complex)` instead
returns a typed `Complex` whose imaginary component is zero.

General Clifford function results must be representable in the selected real
algebra; unsupported or failed evaluation returns a contextual NaN. `abs`,
`signum`, and componentwise rounding on a general multivector without central
scalar meaning return NaN. Scalar-only operations do not discard generator
context.

The `Numeric` supertraits require clone/format/hash, partial comparisons,
arithmetic, and negation. They do not promise `Copy`, `Eq`, `Ord`, or assignment
operators.

### Deferred Numeric methods

These commented signatures are plans, not callable `Number` APIs:

`erf`, `erfc`, `gamma`, `lgamma`, `digamma`, `trigamma`, `tetragamma`,
`elliptic_k`, `elliptic_e`, `zeta`, `exp_polar`, `besselj`, `bessely`,
`besseli`, `besselk`, `polygamma`, `beta`, `zeta_deriv`, `lambertw`,
`hermite`, `assoc_legendre`, `spherical_harmonic`.

## Real

`Real` is a non-exhaustive enum with public `Int(IntType)`,
`Rational(RationalType)`, and `Float(FloatType)` variants. Consumers can construct
variants directly; external matches require a catch-all arm. Use the
normalizing constructors when denominator-one rationals should become integers.

| Methods | Return and behavior |
|---|---|
| `from_int(IntType)`, `from_float(FloatType)`, `from_rational(RationalType)` | `Real`; the rational constructor normalizes denominator one |
| `zero`, `one`, `neg_one`, `nan`, `pi`, `epsilon` | Associated constructors returning `Real` |
| `is_float()` | `bool`; float storage test |
| `is_zero`, `is_one`, `is_neg_one`, `is_negative`, `is_positive` | `bool`; numeric identity/sign tests |
| `is_nan`, `is_infinite`, `is_finite`, `is_int` | `bool`; mathematical value tests |
| `to_f64()` | `f64`; approximate projection |
| `to_int()` | `IntType`; truncating projection, with saturating float cast |
| `to_float()` | `Real`; approximate representation |
| `trunc()` | `Real`; integer result for exact rational inputs, float result for float inputs |
| `total_cmp(&other)` | `Ordering`; real value order, signed zeros equal, NaNs after numbers and ordered by bits |

`Real` implements `Numeric`, `Clone`, `Debug`, `Display`, `Hash`, `PartialEq`, and
`PartialOrd`, with `From<IntType>`, `From<RationalType>` and `From<FloatType>`
conversions that preserve the representation's normalization rules.
Arithmetic supports owned/owned, owned/borrowed, and
borrowed/borrowed operands, plus owned or borrowed unary negation.
`+=` and `-=` accept borrowed `Real` operands and reuse integer or float
accumulators when that representation remains active.
`Real` does not implement `Eq`, `Ord`, or `Copy`.

## RationalType, IntType, and FloatType

`IntType` and `FloatType` are primitive aliases and expose the ordinary `i64` and
`f64` APIs. Primitive integer operators retain Rust's normal primitive behavior;
wrapping integer arithmetic described above belongs to `Number`/`Real` dispatch.

`RationalType` is the public native rational alias. Its fields are private and
always reduced, with a positive denominator.

| Method | Return and behavior |
|---|---|
| `new(numer: i64, denom: i64)` | `RationalType`; either denominator sign; panics for zero or an unrepresentable reduced result |
| `from_parts(numer: i64, denom: u64)` | `RationalType`; full positive unsigned denominator range; panics for zero |
| `from_integer(numer: i64)` | `RationalType`; denominator one |
| `numer()`, `denom()` | `&i64`, `&u64` |
| `is_int()` | `bool`; denominator one |
| `to_integer()` | `i64`; quotient truncated toward zero |
| `trunc`, `fract`, `floor`, `ceil`, `round` | `RationalType`; exact operations; rounding ties away from zero |

The getters borrow the reduced components; `*q.numer()` and `*q.denom()` copy
native values. The references provide read-only access without copying storage.

```rust
use num_anafis::{Number, RationalType, n, r};

let q = RationalType::from_integer(5); // typed rational 5/1
assert_eq!(*q.numer(), 5);
assert_eq!(*q.denom(), 1);
let reciprocal = RationalType::from_parts(1, 5); // typed rational 1/5
assert_eq!(*reciprocal.numer(), 1);
assert_eq!(*reciprocal.denom(), 5);
assert_eq!(Number::from_rational(q), n(5)); // Number normalizes 5/1 to integer
assert_eq!(r(1, 5), n(1) / n(5));
```

`RationalType` implements `Clone`, `Debug`, `Display`, `Hash`, `PartialEq`, `Eq`,
`PartialOrd`, and `Ord`. `+`, `-`, `*`, `/` support owned/owned and
borrowed/borrowed operands; unary negation supports both. It does not implement
`Numeric` or `Copy`. Reduced capacity failures, division by rational zero, and
unrepresentable negation panic. Use `Number` for promotion and elementary
functions; rational operations alone do not promote to floats.

## Complex

`Complex::new(re, im) -> Result<Complex, NumAnafisError>` accepts two
`Into<Number>` arguments and rejects complex/Clifford components. Its public
`re: Real` and `im: Real` fields can be read and changed. The struct is
non-exhaustive, so construct it through its methods, not an external struct
literal. A typed `Complex` remains complex-shaped even with zero imaginary part;
converting it into `Number` applies normalization.

| Inherent methods | Return |
|---|---|
| `zero()`, `one()`, `i()` | `Complex` |
| `conj()` | `Complex` |
| `norm_sq()`, `abs()`, `arg()` | `Real` |
| `add(&rhs)`, `sub(&rhs)`, `mul(&rhs)`, `div(&rhs)` | `Complex` |
| `exp`, `ln`, `sqrt`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh` | `Complex` |

All `Numeric` methods are available through the trait as well. Typed complex
constructors preserve signed zeros for branch evaluation. Arithmetic supports
all owned/borrowed combinations for `+`, `-`, `*`, `/`, and unary negation.
`Complex` implements `Clone`, `Debug`, `Display`, `Hash`, `PartialEq`, and
`PartialOrd`; strict order requires zero imaginary components. Equal complex
values compare `Equal`; distinct non-real values are unordered.

```rust
use num_anafis::{Complex, Number, Real};

let mut z = Complex::new(2, 3).expect("real components");
z.im = Real::from_int(4);
let value: Number = z.into();
assert_eq!(value, num_anafis::c(2, 4));
```

## F32Ext and F64Ext

`F32Ext` is implemented for `f32`; `F64Ext` for `f64`. Methods consume the
primitive and return that same primitive type. Below, `F` means the corresponding
float type, and `Order` means `i32` for `F32Ext`, `i64` for `F64Ext`.
These are real primitive APIs, with no automatic complex or `Number` return.
In `std`, inherent primitive methods take precedence when names overlap;
fully qualified calls such as `F64Ext::sin(x)` select the trait explicitly.

| Methods/signature | Meaning |
|---|---|
| `sin`, `cos`, `tan`, `asin`, `acos`, `atan` | Real trigonometric/inverse functions |
| `atan2(self, other: F)` | Quadrant-aware arctangent |
| `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh` | Real hyperbolic/inverse functions |
| `exp`, `exp_m1`, `ln`, `ln_1p`, `sqrt`, `cbrt` | Elementary functions with primitive naming |
| `powf(self, exponent: F)` | Real float power |
| `floor`, `ceil`, `round`, `trunc`, `fract` | Real rounding/fractional operations |
| `gamma`, `lgamma` | Gamma and `ln(abs(Gamma(x)))` |
| `digamma`, `trigamma`, `tetragamma` | First three logarithmic Gamma derivatives |
| `polygamma(self, order: Order)` | Higher polygamma |
| `erf`, `erfc` | Error function and complement |
| `zeta`, `zeta_deriv(self, order: Order)` | Riemann zeta and its derivative |
| `bessel_j`, `bessel_y`, `bessel_i`, `bessel_k` with `(self, order: Order)` | Integer-order Bessel functions |
| `beta(self, other: F)` | Beta function |
| `elliptic_k`, `elliptic_e` | Complete elliptic integrals using parameter **m = k²** |
| `hermite(self, degree: Order)` | Physicists' Hermite polynomial |
| `lambert_w0`, `lambert_wm1` | Principal and secondary real Lambert W branches |
| `assoc_legendre(self, degree: Order, order: Order)` | Associated Legendre function, with `self` as its argument |
| `spherical_harmonic(self, degree: Order, order: Order, phi: F)` | Current real-valued angular kernel; `self` is theta |

The spherical-harmonic implementation uses a cosine angular factor for both
signs of the order, with its documented phase. It does not return the usual
complex spherical harmonic or a separate sine member for negative order.
Elliptic arguments are parameters, despite the current extension-trait comments
using modulus-style `K(k)`/`E(k)` notation.

```rust
use num_anafis::F64Ext;

let gamma = F64Ext::gamma(5.0_f64);
assert!((gamma - 24.0).abs() < 1e-12);
```

## GeneratorSet (clifford)

Defines generator labels and their diagonal metrics. Generator positions refer
to the sorted label order. Metrics must be `+1`, `-1`, or `0`. Library IDs 0–255
are reserved; custom generators should start at 256.

| Method | Return and behavior |
|---|---|
| `empty()` | `GeneratorSet`; ordinary scalar algebra |
| `new(Vec<(u32, i8)>)` | `Result<GeneratorSet, NumAnafisError>`; sorts labels, rejects duplicates, invalid metrics and unrepresentable blade counts |
| `len()`, `is_empty()` | `usize`, `bool` |
| `id_at(position)`, `metric_at(position)` | `u32`, `i8`; panics outside the generator range |
| `complex_blade()` | `Option<usize>`; oriented full pseudoscalar mask when central and square −1 |
| `union_with(&other)` | `(GeneratorSet, Vec<u8>, Vec<u8>)`; union plus old-position-to-new-position maps |

`union_with` panics for conflicting metrics on a shared label or a union beyond
the platform blade-count capacity (`n >= usize::BITS`). `GeneratorSet` implements `Clone`, `Debug`, `PartialEq`,
`Eq`, and `Hash`. Its fields are private.

## CliffordNumber (clifford)

Dense multivector with real blade coefficients. Blade indices are bitmasks of
generator positions; blade zero is the scalar and blade count is `2^n`.

| Constructor | Return and behavior |
|---|---|
| `zero(gens)` | `Result<CliffordNumber, NumAnafisError>` |
| `scalar(gens, value: impl Into<Number>)` | Checked real/complex scalar embedding; complex `i` maps to the declared central unit |
| `generator(gens, index: u8)` | Checked basis vector at a generator position |
| `from_coeffs(gens, impl IntoIterator<Item = Real>)` | Exactly `2^n` active coefficients in blade-mask order; accepts owned arrays, vectors and iterators |
| `from_sparse(gens, impl IntoIterator<Item = (usize, Real)>)` | Blade-mask terms; repeated masks are added, unspecified blades are zero |
| `nan(gens)` | `CliffordNumber`; panics if the generator count is unrepresentable |

All constructors except `nan` return `Result`. Complex scalar embedding can fail
with `UnsupportedComplexStructure`; blade coefficients themselves never accept
complex values.

| Method | Return and behavior |
|---|---|
| `generator_set()` | `&GeneratorSet` |
| `n_generators()`, `blade_count()` | `usize` |
| `coeff(blade)` | `&Real`; panics outside active blades |
| `set_coeff(blade, value: impl Into<Real>)` | Mutates a real coefficient; panics for an invalid index |
| `coeffs_slice()`, `coeffs_mut_slice()` | `&[Real]`, `&mut [Real]`; active blades only |
| `nonzero_blades()` | Iterator of `(usize, &Real)` |
| `geometric_mul(&other)` | `CliffordNumber`; same operation as `*` |
| `outer_product(&other)` | `CliffordNumber`; wedge product |
| `left_contraction(&other)` | `CliffordNumber`; For homogeneous grades `r <= s`, selects grade `s-r` of the geometric product; zero for `r > s` |
| `scalar_product(&other)` | `Real`; scalar part of the geometric product |
| `grade(k: u32)` | `CliffordNumber`; grade projection |
| `reverse()`, `grade_involution()`, `clifford_conjugate()` | `CliffordNumber`; geometric involutions |
| `norm_sq()` | `Real`; scalar part of `A * clifford_conjugate(A)`, generally signed |
| `geometric_inverse()` | `CliffordNumber`; returns contextual NaN on failed inversion |
| `is_zero`, `is_one`, `is_neg_one`, `is_int`, `is_negative`, `is_positive`, `is_finite` | `bool`; nontrivial integrality/sign require a real grade-zero scalar |
| `is_int_ring()` | `bool`; every active coefficient is integer-valued, including integral floats |
| `to_float()` | `CliffordNumber`; approximate every coefficient |
| `approx_eq_number(&other, &tolerance)` | `bool`; finite maximum coefficient distance, matching algebras |
| `total_cmp(&other)` | `Ordering`; generator context then coefficients |
| `max(&other)`, `min(&other)` | `CliffordNumber`; real scalar order; non-real/unordered inputs return NaN |

`Numeric` is implemented. Operations preserve generator context. Compatible
complex scalar meaning is recognized by `Number`; general unsupported scalar
functions return NaN. `Number::conj` may select a different ambient involution
from the explicitly named `clifford_conjugate` to send the declared imaginary
unit to its negative.

Operators `+`, `-`, `*`, `/`, and unary `-` support owned/borrowed multivectors.
Assignment operators `+=`, `-=`, `*=`, `/=` accept owned or borrowed multivectors.
Division is **right division**: `A / B = A * inverse(B)`. Different generator
sets are merged for products; incompatible shared metrics panic. Multiplication
by a real `i64`/`f64` on the right scales coefficients. Typed `Real` supports
addition, subtraction, and multiplication on either side, and borrowed division
`&mv / &real` or `&real / &mv`.

`CliffordNumber` implements `Clone`, `Debug`, `Display`, `Hash`, `PartialEq`, and
`PartialOrd`. Direct equality includes generator context; `Number` additionally
recognizes central scalar meaning across representations. Direct partial
strict ordering requires matching algebras and real grade-zero scalar operands.
Equal multivectors compare `Equal`; distinct non-real values and NaNs are
unordered. `total_cmp` supplies structural sorting.

```rust
use num_anafis::{Cga, CliffordNumber, Number, Numeric, Real, i};

let mut vector = Cga::e1();
vector.set_coeff(1, 2_i64);
assert_eq!(vector.coeff(1), &Real::Int(2));
let scalar = CliffordNumber::scalar(Cga::gens(), i()).expect("central CGA unit");
assert_eq!(scalar, Cga::ci());
let literal: Number = vector.into();
assert!(literal.is_clifford());
let nilpotent = Cga::eps();
let one = CliffordNumber::scalar(Cga::gens(), 1_i64).expect("real scalar");
assert_eq!(nilpotent.exp(), one + nilpotent);
```

## Cga (clifford)

All unit constructors return `CliffordNumber` in the same Cl(4,1) algebra.
`gens()` returns its `GeneratorSet`, ordered `e1,e2,e3,e+,e−`.

| Constructor | Definition/meaning |
|---|---|
| `e1`, `e2`, `e3` | Euclidean basis vectors, square +1 |
| `e_plus`, `e_minus` | Conformal extra vectors, square +1 and −1 |
| `orig` | `(e_minus - e_plus)/2`, null origin |
| `inf` | `e_minus + e_plus`, null infinity |
| `eps` | Same element as `inf`, dual nilpotent unit |
| `qi`, `qj`, `qk` | `−e2e3`, `−e3e1`, `−e1e2`; Hamilton quaternion orientation |
| `pseudo3d` | `e1e2e3`, Euclidean volume, square −1; not central in full CGA |
| `pseudo5d` | Full oriented pseudoscalar, central and square −1 |
| `ci` | Same element as `pseudo5d`, ordinary complex embedding |
| `sj` | Same element as `e_plus`, split-complex unit, square +1 |

Quaternion units obey `qi*qj=qk` and cyclic variants. The norms of `a+b*ci`,
quaternions, `a+b*sj`, and `a+b*eps` are respectively `a²+b²`, a sum of four
squares, `a²−b²`, and `a²`. Catalog units interact through their actual geometric
products; different labels are not interchangeable imaginary conventions.
`Cga` is a `Clone + Copy + Debug` namespace with no public fields.

## FastClifford<P,Q,R> (clifford)

`P`, `Q`, and `R` are counts of positive, negative, and null generator metrics.
`P+Q+R` must be at most five; using unsupported signatures triggers the size
assertion. Fields are private.

| Method/conversion | Return and behavior |
|---|---|
| `zero()` | `FastClifford<P,Q,R>` |
| `active_generators()` | `usize`; associated signature count |
| `coeffs_slice()`, `coeffs_mut_slice()` | Active `&[Real]` / `&mut [Real]` |
| `total_cmp(&other)` | `Ordering`; coefficient order |
| `FastClifford::try_from(&dynamic)` | `Result<Self, NumAnafisError>`; checks count and metric order |
| `CliffordNumber::from(&fixed)` | Dynamic representation with position-based labels 0..n |

Conversion from a dynamic value keeps ordered coefficients; it checks signature
rather than preserving arbitrary generator labels. Going back to dynamic uses
new position labels. Do not assume a round trip preserves a labeled CGA context.
To make a `Number`, convert through `CliffordNumber` first.

`Numeric`, owned/borrowed `+`, `-`, `*`, `/`, unary negation, and assignment
operators are implemented. Products and multiplication assignment use the fixed
signature kernel; division and elementary functions use the dynamic evaluator.
`Clone`, `Debug`, `Display`, `Hash`, `PartialEq`, and `PartialOrd` are implemented;
strict numerical ordering requires real scalars. Equal values compare `Equal`;
distinct non-real values and NaNs are unordered. `Copy`, `Eq`, and `Ord` are not
implemented.

```rust
use num_anafis::{CliffordNumber, FastClifford, Number, Real};

let mut fixed = FastClifford::<1, 0, 0>::zero();
*fixed.coeffs_mut_slice().get_mut(1).expect("first generator") = Real::Int(2);
let dynamic = CliffordNumber::from(&fixed);
let recovered = FastClifford::<1, 0, 0>::try_from(&dynamic).expect("matching signature");
assert_eq!(fixed, recovered);
let literal: Number = dynamic.into();
assert!(literal.is_clifford());
```

## NumAnafisError

Non-exhaustive, `Clone + Debug + PartialEq + Eq + Display`. With `std`, it also
implements `std::error::Error`. External matches require a catch-all arm.

| Variant | Meaning/details |
|---|---|
| `ExpectedReal` | Typed real component required |
| `ExpectedScalar` | Ordinary scalar meaning required |
| `UnsupportedComplexStructure` | No declared central square-minus-one unit |
| `ActiveGeneratorsExceedSignature(details)` | `details.active`, `details.available` |
| `GeneratorIndexOutOfRange(details)` | `details.index`, `details.active` |
| `DuplicateGenerator { id }` | Repeated generator label |
| `InvalidGeneratorMetric { id, metric }` | Diagonal metric outside `{-1, 0, 1}` |
| `DenseCoefficientLengthMismatch(details)` | `details.expected`, `details.found` |
| `BladeIndexOutOfRange { blade, count }` | Sparse blade mask outside `0..count` |
| `ActiveGeneratorsTooLargeForPlatform { active }` | Blade indexing exceeds platform capacity |
| `MismatchedGeneratorSet` | Incompatible fixed-signature conversion |

The boxed detail structs are not separately re-exported by `lib.rs`; their public
fields remain readable from matched variants. Unsigned `Number::try_from` uses
the standard `TryFromIntError`, not `NumAnafisError`. Rational capacity failures
and invalid primitive projections panic rather than returning this error enum.

## What lib.rs does not expose

There is no public `Scalar`, `Coefficient`, `NumberRepr`, `IntoNumber`, backend
math facade, matrix/eigensolver type, spectral AST/Jet interpreter, or raw storage
enum. Access the representations through the documented types and methods above.
Python classes have their own binding surface, separate from these Rust exports.

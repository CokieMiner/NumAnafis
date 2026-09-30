# NumAnafis

Numeric values for the SymbAnaFis CAS and evaluator. One `Number` holds integers,
rationals, floats, complex values, and optional Clifford multivectors.

## Installation

```toml
[dependencies]
num_anafis = { git = "https://github.com/CokieMiner/NumAnafis" }
```

| Feature | Usage |
|---|---|
| `std` | Enabled by default; disable for `no_std` with an allocator |
| `clifford` | Enable Clifford values, geometric products, and the CGA catalog |
| `python` | Enable Python bindings; requires `std` |

Add `features = ["clifford"]` to the dependency to use geometric algebra.

## Numbers and functions

Import `Numeric` to call elementary functions.

```rust
use num_anafis::{Number, Numeric, c, i, n, r};

let value: Number = n(2) / n(3);
assert_eq!(value, r(2, 3));
assert_eq!(r(4, 9).sqrt(), r(2, 3));
assert_eq!(r(2, 3).pow(&n(-3)), r(27, 8));
assert_eq!(r(-4, 9).sqrt(), c(0, r(2, 3)));
assert_eq!(i() * i(), n(-1));
assert_eq!(5 * i(), c(0, 5));
assert!(n(2.0).is_float());
assert!(n(2.0).is_int());
```

| Constructor | Result |
|---|---|
| `n(value)` | A `Number` from a primitive or numeric value |
| `r(numer, denom)` | Exact rational; a zero denominator returns NaN |
| `c(re, im)` | `re + i*im` |
| `i()` | Imaginary unit |

Standard `From`/`Into` conversions are also supported: `let x: Number = 2.into();`.
Use `Number::try_from(value)` for `u64`/`usize`; out-of-range inputs return an error.

`n(2)` creates an integer; `n(2.0)` keeps the float representation. Exact division,
integer powers, and perfect roots stay exact when representable. Rationals with
denominator one reduce to integers. Other function results may be approximate.
For typed components, use `Real`, `Complex::new`, and
`RationalType::from_parts(numer, denom)`; the latter accepts a positive `u64`
denominator.

`Numeric` supports trigonometric, hyperbolic, exponential, logarithmic, root,
power, rounding, sign, and unnormalized `sinc` functions. Complex logarithms and
noninteger powers use principal branches; `cbrt` keeps real cube roots for real
inputs. Special functions such as Gamma are currently available on primitive
`f32`/`f64` through `F32Ext`/`F64Ext`, rather than on `Number`.

## Inspection and conversion

- `is_int()` tests whether the represented value is a real integer. Finite
  integral floats pass; the predicate does not establish computation exactness.
- `is_int_ring()` checks that every real component or blade coefficient is an integer.
- `re()` and `im()` project the scalar and declared complex-unit coefficients.
  In CGA, the imaginary unit is I5.
- `conj()`, `norm_sq()`, and `arg()` provide conjugation, squared norm, and argument.
  General multivectors have a signed norm and no scalar argument.
- `to_f64()` and `to_int()` project the real scalar component; `to_int()` truncates.
  `to_float()` converts every component to the float representation.
- `extract_complex()` converts the entire value to typed complex components,
  rejecting additional blades. Primitive projections require a scalar value.
- `max()` and `min()` compare real scalars; non-real inputs return NaN.
  `total_cmp()` supplies an explicit sorting convention.

Finite equality compares represented values across numeric tiers. NaN remains
unordered and unequal to itself. `Number`, `Real`, `Complex`, and Clifford
values do not implement `Eq`.

The native backend uses `i64`, `Ratio<i64,u64>`, and `f64`. Integer arithmetic
wraps on overflow. Rational arithmetic panics when a reduced result exceeds
native capacity. No operation automatically increases backend capacity.

## Clifford values

With the `clifford` feature:

```rust
use num_anafis::{Cga, Number, Numeric, i, n};

let imaginary: Number = Cga::ci().into();
assert_eq!(imaginary, i());
assert_eq!(&imaginary * &imaginary, n(-1));
let mixed: Number = Cga::e1().into();
assert!((mixed + i()).is_clifford());

let translation = Cga::eps();
let translated = translation.exp();
assert_eq!(translated.coeff(0), &num_anafis::Real::one());
```

Blade coefficients are real. Complex scalars map to a central pseudoscalar whose
square is −1; checked constructors reject algebras without that declared unit.
An unrepresentable function result returns NaN. Ordinary operations preserve
Clifford generator context.

In CGA, `ci()` is `pseudo5d()`. `pseudo3d()` is the Euclidean volume element.
`CliffordNumber::norm_sq()` is a signed geometric form for general multivectors;
use coefficient distance for approximate comparison. `Number::conj()` preserves
ordinary complex conjugation on its declared scalar subalgebra; the named
`clifford_conjugate()` retains its geometric definition.

See [Clifford usage and unit definitions](src/clifford/README.md) for the catalog,
products, and norms.

See [the complete public Rust API guide](API.md) for all exports and CAS integration.

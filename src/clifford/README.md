# Clifford algebras

Enable the `clifford` feature to use `CliffordNumber`, `GeneratorSet`, `Cga`, and
`FastClifford<P,Q,R>`. Blade coefficients are real integers, rationals, or floats.
A multivector can be converted into `Number` for use in the evaluator.

```rust
use num_anafis::{Cga, Number, Numeric, Real};

let bivector = Cga::qi() * 0.3_f64;
let rotor = bivector.exp();
let value: Number = rotor.into();

let mut vector = Cga::e1();
vector.set_coeff(1, 2_i64);
assert_eq!(vector.coeff(1), &Real::Int(2));
```

`*` is the geometric product. `outer_product` computes the wedge product;
`left_contraction`, `scalar_product`, and `grade` provide the corresponding
projections. `coeffs_slice` and `coeffs_mut_slice` expose active coefficients.
`GeneratorSet::new` sorts labels and validates metrics. `from_coeffs` accepts
exactly `2^n` real coefficients; `from_sparse` accepts blade-mask terms and adds
repeated masks. Both return `Result` for invalid input.
Numerical comparison and `min`/`max` require real scalar operands; extrema return
NaN for non-real or unordered values. `total_cmp` supplies deterministic sorting.
`FastClifford<P,Q,R>` specifies counts of positive, negative, and null generators
and supports up to five generators. Its algebraic products avoid dynamic algebra
construction; elementary functions and division use the dynamic evaluator.

## CGA catalog

The ordered signature is `(e1²,e2²,e3²,e+²,e−²) = (+1,+1,+1,+1,−1)`.

| Constructor | Definition | Meaning |
|---|---|---|
| `pseudo3d()` | `e1 e2 e3` | Oriented Euclidean volume; square −1 |
| `pseudo5d()`, `ci()` | `e1 e2 e3 e+ e−` | Central complex unit; square −1 |
| `qi()`, `qj()`, `qk()` | `−e2 e3`, `−e3 e1`, `−e1 e2` | Hamilton quaternion units |
| `sj()` | `e+` | Split-complex unit; square +1 |
| `orig()` | `(e− − e+)/2` | Null origin |
| `inf()`, `eps()` | `e− + e+` | Null infinity and dual unit |

`ci()` aliases I5 so it commutes with every CGA blade. I3 commutes with the
Euclidean Cl3 subalgebra but anticommutes with `e+` and `e−`; it remains a distinct
geometric element. Quaternion orientation gives `qi*qj=qk`, `qj*qk=qi`, and
`qk*qi=qj`. `orig·inf=−1`.

The spans of `1,sj` and `1,eps` obey split-complex and dual-number laws. Their
units commute with Euclidean quaternion bivectors and anticommute with Euclidean
vectors. Mixing catalog units follows this embedding's blade products.

`reverse()` is geometric reversion. `clifford_conjugate()` is reversion followed
by grade involution. `norm_sq()` returns the scalar part of
`A*clifford_conjugate(A)`: central CGA complex/quaternion values have their usual
squared modulus, split-complex values have `a²−b²`, and dual values have `a²`.
For general multivectors this is a signed quadratic form. General conjugate
products need not be scalar or provide an inverse. Approximate comparison uses
maximum real coefficient distance.

## Complex promotion and functional calculus

Blade coefficients are always real. `CliffordNumber::scalar` maps `a+bi` to
`a+bJ`, where J is the oriented full pseudoscalar only if it is central and
squares to −1. Checked construction rejects complex values in other algebras.
`Number::re()` projects the grade-zero coefficient. `Number::im()` projects the
coefficient of J, or returns zero if the algebra has no declared complex unit.
These projections exclude other blades. `Number::extract_complex()` rejects
additional blades and returns a typed complex value without generator context.
Compatible scalar values are remapped when entering another algebra. Clifford
arithmetic and functional calculus preserve generator context.

Elementary functions return NaN when the result cannot be represented in the
algebra. Nilpotent exponentials with `B²=0` preserve the exact result `1+B`.

For a declared central complex structure, `Number::conj` uses an anti-involution
that sends its imaginary unit to its negative. In ranks congruent to 3 modulo 4
this is reversion; in ranks congruent to 1 it is Clifford conjugation. The named
`CliffordNumber::clifford_conjugate` always retains its geometric definition.
`Number::norm_sq` uses the corresponding scalar conjugate product, so scalar
complex meaning and the ambient involution remain consistent.

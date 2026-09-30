# Special-function accuracy work

Evidence: precision run `5ff1fa1381fe`, dated 2026-09-30. The run evaluated
47 primitive methods at both widths, with 1,000 random samples, 100 targeted
samples, and represented edge cases per function/width; seed 42.
Native and reference evaluations used two-second work budgets.

The [measurement report](../../tools/results/analyze.md) records the results.
The [precision CLI](../../tools/precision_suite/README.md) describes domains,
reference conventions, sampling, and error calculations. SQLite data is local
and ignored by Git. These observations describe this run, not subsequent edits.

ULP error uses spacing at the unrounded reference. Agreement between two
reference precisions is empirical, not certified rounding. Large ULP errors
near roots remain failures of ULP accuracy even when absolute error is small.
Reference errors and timeouts are unscored; they do not establish a kernel bug.

## 1. Domain, limiting-value, and termination failures

Resolve these before optimizing approximations. Reproduce each case separately
and verify its mathematical domain, target rounding, and signed-zero convention.

| Implementation | Observed issue | Work required |
|---|---|---|
| [lambert_w.rs](lambert_w.rs) | Both branches produce NaNs for valid inputs; near the branch point some results collapse to -1. | Review initial estimates, branch-point handling, and convergence acceptance. |
| [bessel_jy.rs](bessel_jy.rs), [bessel_ik.rs](bessel_ik.rs) | J and I return zero for representable nonzero small-argument results. | Review small-argument evaluation and recurrence scaling; distinguish genuine underflow from an unstable calculation. |
| [polygamma.rs](polygamma.rs) | Tetragamma returns a finite constant at minimum subnormal inputs; digamma and general polygamma also have exceptional-value mismatches. | Review near-zero singular terms, reflection, and intermediate overflow. |
| [gamma.rs](gamma.rs) | Maximum finite positive inputs exceed the work budget; positive infinity produces NaN. | Ensure bounded large-argument evaluation and handle defined infinite limits. |
| [bessel_ik.rs](bessel_ik.rs) | K produces NaN where the mathematical value rounds to positive infinity. | Review small-positive-argument scaling and overflow handling. |
| [hermite.rs](hermite.rs) | Extreme finite inputs and infinities can produce NaN instead of the expected signed infinity. | Review polynomial limiting values and overflow through the compensated recurrence. |

Concrete measured f64 cases:

- `lambert_w0(-5e-324)` returned `2.465190328815662e-32`;
  the reference rounds to `-5e-324`.
- `lambert_wm1(-5e-324)` returned NaN; the reference is approximately
  `-751.0615595398791`.
- `bessel_j(-4.6216115207296995e-32, 8)` returned zero; the reference is
  approximately `2.0164479138637223e-258`, a representable nonzero result.
- `bessel_i(-6.291736556510495e-30, -2)` returned zero; the reference is
  approximately `4.9482436120663175e-60`.
- `tetragamma(5e-324)` returned approximately `-2.4041138063191885`;
  the mathematical value rounds to negative infinity.
- `bessel_k(5e-324, 7)` returned NaN; the mathematical value rounds to
  positive infinity.

Beta pole signs, signed zeros, and other exceptional-value mismatches require
explicit contracts before being classified as implementation defects.
In particular, a reference failure at negative infinity is not evidence that
the native result is wrong.

## 2. Accuracy on ordinary sampled inputs

Percentages below use only finite, scored random results, excluding targeted
and edge samples. They are sample frequencies, not domain-wide failure rates.
Gamma has fewer finite results because many sampled positive inputs overflow
the target format. Worst errors in this table use the same random subset.

| Function | f32 results above 10 ULP | f64 results above 10 ULP | Worst random f64 error |
|---|---:|---:|---:|
| Gamma | 51% | 75% | 515 ULP |
| Bessel Y | 47% | 52% | 8,875 ULP |
| Bessel J | 42% | 39% | 231,907 ULP |
| Beta | 40% | 34% | 50 ULP |
| Zeta derivatives | 16% | 32% | 1,947 ULP |
| erfc | 13% | 15% | 51 ULP |
| Spherical harmonics | 10% | 10% | 5,990 ULP |

Investigate the following numerical regimes:

- **Gamma:** recurrence accumulation for positive arguments and reflection
  for negative arguments. Keep accuracy work separate from termination fixes.
- **Bessel J/Y:** recurrence direction, transition between approximation
  regimes, roots, and large-argument phase reduction.
- **Beta:** logarithmic evaluation, reflection, cancellation, and the sign of
  the result for negative arguments.
- **Zeta derivatives:** differentiation order, quadrature conditioning,
  distance to poles, and agreement with independent derivative references.
- **erfc:** positive tails, especially arguments around 8–10; small absolute
  errors there do not imply small ULP errors.
- **Spherical harmonics:** propagation through associated Legendre values,
  normalization, and the trigonometric factor near its zeros. Preserve the
  declared cosine projection and negative-order phase convention.

## 3. Smaller accuracy issues and stronger results

| Function | Approximate worst finite error across f32/f64 | Remaining attention |
|---|---:|---|
| erf | 2.35 ULP | No exceptional-value mismatches recorded; improve approximation accuracy if the target is below this. |
| Elliptic K | 2.38 ULP | No exceptional-value mismatches recorded; its endpoint pole was unscored. |
| Trigamma | 3.55 ULP | One reference failure per width; verify that case separately. |
| Bessel K | 4.58 ULP | Finite accuracy is relatively strong; exceptional-value failures remain urgent. |
| Hermite | 0.5 ULP | Finite samples performed well; do not overlook the limiting-value failures above. |
| Elliptic E | 14.93 ULP | Review difficult parameter regimes and the negative-infinity limit. |
| Associated Legendre | 98.19 ULP | Review recurrence near roots/endpoints and signed zeros. |
| lgamma | 270.24 ULP | Review negative-argument roots and reflection; verify exceptional-value policy. |
| Digamma | 397.59 ULP | Review roots and small negative arguments. |
| General polygamma | 15,854 ULP | Review order-dependent reflection and singular behavior. |
| Tetragamma | 15,794,528 ULP | Review singular behavior first; finite f32 errors also need attention. |
| Zeta | 56.33 ULP | Review negative arguments, roots, and reflection. |

No special function matched the rounded reference in every scored case at
both widths. The primitive run does not independently validate complex helper
kernels or deferred complex-domain special-function integration.

## 4. Related primitive finding

`atanh(-0.9999999999999999)` returned approximately `-18.36840028483855`
instead of `-18.714973875118524`: absolute error approximately `0.3466`.
This is outside `special/`; isolate the primitive extension/backend path in
[f64ext.rs](../ext/f64ext.rs) and check f32 separately before changing a kernel.

## Verification of repairs

Preserve raw input bits and integer parameters in regressions. Check numerical
error, target-rounded exceptional values, and signed zero separately. Include
adjacent inputs at roots, poles, branch points, and approximation boundaries.
Do not infer accuracy from a fixed absolute EPSILON threshold or from a plot.
Record accuracy targets by regime and evaluate std and no_std paths separately.

```sh
python3 -m tools.precision_suite run --function gamma,bessel_j,bessel_y --backend both --samples 1000 --targeted 100 --seed 42 --timeout 2 --reference-timeout 2
python3 -m tools.precision_suite report --only-severe
```

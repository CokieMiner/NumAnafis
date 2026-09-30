# MP backend integration

MP provides UInt, Int, Rational, Float, and complex values over the three real
tiers. Number remains the public evaluator value; Real holds Int, Rational, or
Float. Complex components and Clifford blade coefficients remain real-only.
UInt supplies unsigned magnitudes and positive rational denominators.

## Backend work

- Adapt `src/types/int.rs`, `rational.rs`, and `float.rs` to the MP types. Use
  existing backend gcd, quotient/remainder, and root/remainder operations;
  cross-cancel rational operands before multiplication.
- Replace fixed-width intermediates and Copy assumptions with borrowed
  arithmetic, owned values, and reusable limb/scratch buffers. Adapt native
  const constructors where owned representations require runtime construction.
- Generalize `src/types/real_compare.rs` using backend integer magnitudes and
  float significand/exponent access. Preserve exact cross-tier equality/hash,
  unordered arithmetic NaNs, and branch-significant signed zeros.
- Traverse integer exponent magnitude bits without an i64 projection.
- Construct constants and perform approximate calculations at backend working
  precision. Preserve exact integer/rational coefficients; reserve primitive
  conversions for explicit API boundaries, including Python bindings.
- Reuse MP complex kernels where tier dispatch permits. Preserve mixed real
  components and exact promotion; keep scalar imaginary-unit mapping consistent
  with the declared central Clifford complex blade (I5 in CGA).

## Numerical accuracy

- Document domains, branch conventions, convergence assumptions, and error
  measures. Distinguish proven bounds, numerical estimates, and heuristics.
- Native Taylor margins (64 extra terms, 8 tail terms), QR cap (200 iterations
  per dimension), residual/deflation factors (100), and sqrt(epsilon) clustering
  are heuristics. Reassess them using backend precision, rounding, and
  conditioning; do not transfer them unchanged as accuracy guarantees.
- Use function-specific Taylor remainder bounds where available. Check the
  convergence domain and account separately for rounding, cancellation, and
  guard bits. A small terminal window does not certify the remaining error.
- Preserve exact nilpotent termination: strictly upper triangular n-by-n N
  satisfies N^n=0, so analytic functional calculus needs terms through n-1.
- Derive finite-difference steps from truncation and roundoff analysis. The
  central-difference exponent 1/3 has this justification; its local scale still
  depends on the function.
- Check decomposition/solve residuals and conditioning. Distinguish backward
  error from forward accuracy, and numerical rank from exact algebraic rank.
- Keep work caps separate from accuracy acceptance. Decline evaluation or use a
  justified fallback when work is exhausted or accuracy criteria are unmet.

Taylor remainder reference: [Davies–Higham, section 2 and Algorithm 2.6](https://eprints.maths.manchester.ac.uk/156/1/paper.pdf).

## Validation

- Test native and MP backends, std/no_std, Clifford interoperability, and
  primitive/Python conversion boundaries.
- Compare multiple precisions against exact cases and independent references.
  Cover clustered spectra, defective matrices, large dynamic ranges,
  cancellation, and branch boundaries; record achieved errors and convergence.
- Benchmark scalar arithmetic, Clifford products, closed-form exponentials,
  inversion, and spectral evaluation separately, including allocations.
- Run the repository's Clippy, documentation, structure, and import checks.
  Tests supplement mathematical justification; they do not establish bounds.

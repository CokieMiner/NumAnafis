## Interval

Interval arithmetic `[lo, hi]` for rigorous error bounding.

- **Rounding**: outward rounding on every operation (lo rounds down, hi rounds up)
- **All `Numeric` trait functions**: must be monotonicity-aware
  - Monotone increasing (exp, sinh): `[f(lo), f(hi)]`
  - Monotone decreasing (exp_neg): `[f(hi), f(lo)]`
  - Non-monotone (sin, cos): requires critical point analysis
- **Intersection / union** operations
- **Width / midpoint** queries

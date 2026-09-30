# Primitive accuracy CLI

Compare all `F32Ext` and `F64Ext` methods directly with high-precision mpmath
references. Rust must be installed. Python dependencies are `mpmath`, `numpy`,
and `matplotlib`; reports and exports do not load plotting libraries.
Run commands from the repository root.

```sh
python3 -m tools.precision_suite run --list
python3 -m tools.precision_suite run --function gamma,lgamma,erf --backend both --samples 1000 --seed 42
python3 -m tools.precision_suite report --only-severe
python3 -m tools.precision_suite plot --function lgamma
```

The default database is `tools/results/precision.sqlite3`. Each invocation
creates a separate run. `report`, `plot`, and `export` accept `--db` and `--run`;
otherwise they select the latest run. Plots are PNG files under
`tools/results/charts/32` and `64`, or the directory supplied through `--output`.

`--mode random` samples the configured real domains. `grid` samples evenly
spaced primary inputs; `bits` samples raw bit patterns. Extra parameters use
deterministic sampling in each mode. All sampled runs also test signed zeros,
subnormals, endpoints, infinities, NaNs, and adjacent represented inputs.
Gamma-family poles and both positive `lgamma` roots receive explicit coverage.
`--targeted` controls additional samples near poles, endpoints, and zero.
Counts refer to evaluations; sampled edge cases may repeat.

`--no-std` evaluates the library's libm paths. `--fail-ulp 1` adds an error
threshold to the exit status. A completed run returns 1 for threshold failures,
exceptional-value or signed-zero mismatches, native failures, or reference
failures; configuration/build errors return 2. Known Gamma/psi integer poles, the zeta pole, and the elliptic-k endpoint remain explicitly unscored.
Native and reference time limits are configured through `--timeout` and
`--reference-timeout`, in seconds. Timeouts are measurements of the configured
work budget, not proofs that an algorithm cannot terminate.

## Error definitions

Inputs are rounded to the selected target before reference evaluation. Raw
input/output bits are preserved, including f32 NaN encodings sent to Rust.
The reference receives those represented values, not their decimal labels.

For a finite reference `r`, precision `p`, and minimum normal exponent `emin`:

```text
spacing(r) = 2^(max(emin, floor(log2(abs(r)))) - p + 1)
ULP error  = abs(output - unrounded reference) / spacing(reference)
```

Zero uses minimum subnormal spacing. Binary32 uses `p=24, emin=-126`; binary64
uses `p=53, emin=-1022`. At a power-of-two boundary, spacing belongs to the
reference's binade. Consequently, an adjacent representable value below 1
is 0.5 ULP away by this convention, while its representable-step distance is 1.
`ulp_distance` counts target representable steps to the nearest-even rounded
reference, treating the two zero encodings as one numerical position.
Signed-zero mismatches are recorded separately.

Absolute and relative errors accompany ULP error; relative error is undefined
at an exact zero. No absolute `EPSILON` threshold hides large ULP errors near
roots. High-precision text preserves errors that Python floats cannot represent.
The rounded reference is computed directly in the target format, including
subnormal ties and overflow; f32 does not round through f64 first.

Reference precision starts at 128 bits and doubles until target rounding agrees
at two successive precisions, up to 1024 bits. Configure these with
`--reference-bits` and `--max-reference-bits`. Agreement is an empirical check,
not certified correct rounding or a proven error bound. Nonreal reference
outputs are recorded as failures rather than projected onto the real axis.
Real-domain rejection, known exceptional values, and poles have separate
reference statuses. Elliptic integrals use parameter `m`; spherical harmonics
use the crate's cosine projection with its negative-order phase and no extra
sqrt(2) factor.

## Exhaustive f32 scans

There are 4,294,967,296 f32 bit patterns, including nonfinite values. Exhaustive
mode streams a contiguous raw-bit range; it never creates the full input list.
For multiargument functions only the primary argument is swept; supply fixed
extra arguments in trait order, for example `bessel_j --extras 2`.

```sh
python3 -m tools.precision_suite run --function sin --backend 32 --mode exhaustive --start-bits 0x3f7ff000 --end-bits 0x3f801000
python3 -m tools.precision_suite run --function bessel_j --backend 32 --mode exhaustive --extras 2 --shards 8 --shard 0 --db tools/results/bessel-j-0.sqlite3
```

The start is inclusive and the end exclusive. Omitting endpoints selects the
full bit space. Shards are disjoint contiguous subdivisions; shard indices start
at zero. Use separate databases for parallel jobs. Throughput and ETA reflect
the measured part of the current scan; reference cost varies across domains,
so extrapolation to all f32 values can be very misleading.

SQLite commits counters, retained cases, and the cursor atomically per batch.
Interrupt with Ctrl-C, or use `--max-batches 10` for a controlled checkpoint.
Resume with the printed run ID and identical options through `--resume ID`.
Resumption validates configuration, native binary hash, reference source hash,
Python/mpmath versions, and platform. Completed runs cannot be resumed.

Complete counters and maxima include every evaluated case. Sampled runs retain
all cases; exhaustive scans retain the worst case in each of 4096 raw-bit bins,
plus 20 worst cases and 20 failure examples per function/width. Exhaustive plots
therefore show retained observations, not every float or a continuous-domain
guarantee. `plot --report` prints observed input-bin summaries. Large sampled
plots use a deterministic reservoir controlled by `--max-points`; integer
parameter facets prioritize the highest retained errors through `--max-slices`.
Log plots explicitly clip display values; stored measurements remain intact.
For extreme f64 ranges, explicitly labeled coordinate transforms avoid plotting
overflow; these transforms do not change error measurements.

```sh
python3 -m tools.precision_suite export --output tools/results/export.json
python3 -m unittest discover -s tools/precision_suite/tests
```

JSON is an optional streaming export, not the primary result store. Nonfinite
numbers use tagged objects such as `{"$float":"inf"}`; arbitrary strings and
signed zeros remain distinct.

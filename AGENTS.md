# NumAnafis - Project Guidelines

This document specifies the architecture, invariants, memory safety contracts, performance discipline, module organization, and verification procedures for the `NumAnafis` Rust crate (`num_anafis` in Cargo).

## 1. Priorities & System Architecture

`NumAnafis` provides the numeric foundation for `SymbAnaFis`:
- **`Number`**: One public evaluator value with direct real, complex, and optional Clifford variants.
- **`Real`**: Integer, rational, and approximate tiers (`i64`, `Ratio<i64,u64>`, `f64`).
- **`Complex`**: Two real-only components; never recursive `Number` components.
- **`Numeric`**: Supported elementary evaluator functions. Deferred special signatures stay in its TODO comment, without production stubs.
- **Clifford**: Real-only blade coefficients. Complex scalars map to the declared oriented central pseudoscalar when it squares to -1. The CGA complex unit is I5.
- **Backend boundary**: Arithmetic facades operate on the native integer, rational, and float aliases.

The objective is to avoid unnecessary approximation, not to solve native overflow. CAS AST algebra and commutativity are consumer responsibilities. Retain mathematical `is_int()` for evaluation decisions.

**Design Priorities (strictly ordered):**
1. Mathematical correctness and representation invariants.
2. Memory safety and sound encapsulation.
3. Architecture portability across 32- and 64-bit targets.
4. Optimal execution efficiency (allocation-free promotion paths, in-place evaluation).
5. Ergonomic public APIs conforming to standard Rust idioms.

Correctness precedes performance optimization. Performance never excuses an unproved invariant, and proved invariants must avoid unnecessary allocation or branching.

**Subsystem Boundaries:**
- **`types/`**: Native backend facades; real arithmetic, comparison/hash, typed complex kernels, `Number` dispatch, and constructors.
- **`traits/`**: The `Numeric` evaluator contract.
- **`clifford/`**: Generator sets, storage, products, catalog constructors, faithful embeddings, and spectral evaluation.
- **`ext/`**: Real primitive extension traits (`F32Ext`, `F64Ext`).
- **`special/`**: Existing real special-function kernels, independent of deferred `Numeric` integration.


---

## 2. Correctness, Invariants & Contracts

- **Representation Invariants**:
  - **Representation**: Promote domains when needed while preserving exact components. Normalize denominator-one rationals and complex imaginary positive zero; never relabel approximate floats as exact integers or implicitly discard small imaginary components.
  - **Rationals**: Reduced signed numerator and strictly positive unsigned denominator. Reduce wide intermediates before narrowing; native capacity failures are explicit.
  - **Equality**: Exact equality of represented finite values across real tiers; hashes must agree. Arithmetic NaNs are unordered and nonreflexive, so `Number`, `Real`, `Complex`, and Clifford values must not implement `Eq`.
  - **Clifford Invariants**: Structural fields are encapsulated; mutable slices expose only active real coefficients. Keep scalar imaginary units consistent across promotion and matrix embedding. Geometric conjugate norms may be signed; error checks use positive coefficient distance.
- **Precondition Validation**: Validate preconditions at public and dispatch boundaries. Internal hot-path kernels must not propagate unreachable `Option`, `Result`, or error flags.
- **Arithmetic Modes**: The native integer backend is explicitly modular: `add`/`sub`/`mul`/negation wrap on overflow, identically in debug and release. Overflow detection and bounded arithmetic must be explicit.
- **Production Completeness**: Placeholder stubs (`todo!()`, `unimplemented!()`) and dead branches are forbidden in production code.

---

## 3. Memory Safety & Performance Discipline

- **SAFETY Proofs**: Every `unsafe` block requires an immediately preceding `// SAFETY:` comment proving the obligations applicable to that operation: bounds, initialization, aliasing, alignment, lifetimes, capacity, or target prerequisites.
- **Lint Allowances**: Every `allow` or `expect` attribute requires an explicit, descriptive `reason = "..."` and the narrowest useful scope. Do not suppress a lint instead of correcting invalid code or a misplaced cfg gate.
- **Buffer Reuse**: Commutative operations may swap operands to reuse larger pre-allocated buffers. In-place arithmetic (`+=`, `*=`, etc.) must avoid reallocating when capacity permits.
- **Hot-Path Auditing**: Eliminate loop-invariant checks, redundant normalization, duplicate calculations, and unnecessary allocations.

---

## 4. Module Organization & Visibility

- **`mod.rs` as Module Registry**: `mod.rs` serves exclusively as a structural module registry. It contains zero code (no structs, enums, functions, constants, or traits). All types reside in dedicated child files.
  Item ordering in `mod.rs`:
  1. `//!` module documentation
  2. Inner attributes
  3. `use super::{...}`
  4. `mod ...;`
  5. `pub use ...;`
  6. `#[cfg(test)] mod tests;` (last, only if tests exist)
- **No Sibling Preludes**: Do not import a sibling's private implementation through a parent's private `use`. Use `super::sibling::Thing` directly or an explicit subsystem facade reexport.
- **File Sizing**: Keep production files focused and normally at most 500 lines. Do not split or pad a file solely to meet a line target.
- **Visibility Gates**: Restrict visibility to the minimum necessary scope. `pub(super)` and `pub(in ...)` are forbidden. Use `pub(crate)` where crate access is needed; plain `pub` inside sealed private modules may support parent reexports without exposing an external API.
- **Import Paths & Grouping**: Parent dependencies use `super::Thing`; direct siblings use `super::sibling::Thing`. Cross-subsystem access uses facades (`crate::subsystem::Thing`). Group imports by blank line in order: `core`/`std`, `alloc`, external crates, `crate::`, `super::`, `self::`/child. Follow rustfmt's ordering within groups.

---

## 5. Source-File Structure & Technical Documentation

### Monotonic Call-Graph Ordering (The Data Cascade)
Source files follow downward execution order matching chronological data flow:
1. **Module Header & Declarations**: Module documentation (`//!`), attributes, imports, constants, and type definitions.
2. **Primary Driver / Facade Entry Point**: Inherent method or driver; validates boundaries and initiates execution.
3. **Algorithmic Sequence**: Functions appear in chronological execution order.
4. **Infallible Leaves & Leaf Kernels**: Arithmetic primitives and raw loops placed at the bottom of the file.

### Technical Documentation & Comments
- **Scientific Progression**: Comments must explain the active mathematical identity, precision guarantees, and representation invariants.
- **Rustdoc Standards**: Document API purpose, mathematical domain, panics, errors, and preconditions.
- **Tone & Formality**: Academic, scientific, and concise. State mechanical invariants and formal proofs.

---

## 6. Test Layout & Verification

- **Strict Test Separation**: Production files contain no test logic. Tests reside in dedicated `tests.rs` or `tests/` directories located directly in the module folder, declared only as the final item in `mod.rs`.
- **Coverage Requirements**: Core operations require property-based tests covering edge values, zero handling, and representation promotion.
- **Library Verification**:
  - `cargo check --lib`
  - `cargo clippy --lib`
  - `cargo test --lib`
  - `python3 tools/structure_audit.py`
  - `python3 tools/import_audit.py`
  - `python3 tools/check_allows.py`
  - `python3 -m unittest discover -s tools/audit/tests`

---

## 7. Scope Discipline & Workflow Limits

- **Scope Discipline**: Modifications must directly improve correctness, performance, or guideline compliance. Unrelated refactoring or aesthetic churn is prohibited.
- **Public API Stability**: Stable public APIs must not be renamed or altered without explicit justification.
- **Audit Adherence**: All changes must maintain zero findings on `tools/structure_audit.py` and `tools/import_audit.py`.

use core::{
    cmp::Ordering,
    hash::{Hash, Hasher},
};
use std::collections::hash_map::DefaultHasher;

#[cfg(feature = "clifford")]
use alloc::vec::Vec;

#[cfg(feature = "clifford")]
use pyo3::PyErr;
use pyo3::{
    Bound, Py, PyAny, PyRef, PyResult, Python,
    exceptions::PyValueError,
    prelude::PyModule,
    pyclass, pyfunction, pymethods, pymodule,
    types::{PyAnyMethods, PyInt, PyModuleMethods},
    wrap_pyfunction,
};

#[cfg(feature = "clifford")]
use crate::clifford::{Cga, CliffordNumber, GeneratorSet};
#[cfg(feature = "clifford")]
use crate::error::NumAnafisError;
use crate::{
    traits::Numeric,
    types::{Number, n as raw_n, r as raw_r},
};

#[cfg(feature = "clifford")]
fn map_error(e: &NumAnafisError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

// ============================================================================
// PyNumber
// ============================================================================

#[pyclass(name = "Number", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyNumber(Number);

#[pymethods]
impl PyNumber {
    #[new]
    fn new(x: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(existing) = x.extract::<PyRef<'_, Self>>() {
            return Ok(existing.clone());
        }
        if x.is_instance_of::<PyInt>() {
            // Native capacity failures remain explicit; integer construction
            // must not quietly round through a float.
            x.extract::<i64>().map(|value| Self(raw_n(value)))
        } else {
            x.extract::<f64>().map(|value| Self(raw_n(value)))
        }
    }

    #[staticmethod]
    fn from_ints(num: i64, den: i64) -> Self {
        Self(raw_r(num, den))
    }

    #[staticmethod]
    const fn epsilon() -> Self {
        Self(Number::epsilon())
    }

    fn __repr__(&self) -> String {
        format!("Number({})", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }

    fn __int__(&self, py: Python) -> PyResult<Py<PyAny>> {
        if !self.0.is_int() {
            return Err(PyValueError::new_err(
                "Number is not an integral real value",
            ));
        }
        let builtins = py.import("builtins")?;
        let convert = builtins.getattr("int")?;
        if self.0.re().is_float() {
            // Python converts the represented binary float exactly. Parsing its
            // shortest decimal display could change integers beyond 2^53.
            convert.call1((self.0.to_f64(),)).map(Bound::unbind)
        } else {
            convert.call1((self.0.to_int(),)).map(Bound::unbind)
        }
    }

    fn __float__(&self) -> f64 {
        self.0.to_f64()
    }

    fn __bool__(&self) -> bool {
        !self.0.is_zero()
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.0.hash(&mut hasher);
        hasher.finish()
    }

    fn to_int(&self, py: Python) -> Option<Py<PyAny>> {
        self.__int__(py).ok()
    }

    // --- arithmetic ---

    fn __add__(&self, other: &Self) -> Self {
        Self(&self.0 + &other.0)
    }
    fn __sub__(&self, other: &Self) -> Self {
        Self(&self.0 - &other.0)
    }
    fn __mul__(&self, other: &Self) -> Self {
        Self(&self.0 * &other.0)
    }
    fn __truediv__(&self, other: &Self) -> Self {
        Self(&self.0 / &other.0)
    }
    fn __neg__(&self) -> Self {
        Self(-&self.0)
    }
    fn __abs__(&self) -> Self {
        Self(self.0.abs())
    }
    fn __pow__(&self, exp: &Self, _mod: Option<&Self>) -> Self {
        Self(self.0.pow(&exp.0))
    }

    // --- comparisons ---

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __ne__(&self, other: &Self) -> bool {
        self.0 != other.0
    }
    fn __lt__(&self, other: &Self) -> bool {
        self.0.partial_cmp(&other.0) == Some(Ordering::Less)
    }
    fn __le__(&self, other: &Self) -> bool {
        matches!(
            self.0.partial_cmp(&other.0),
            Some(Ordering::Less | Ordering::Equal)
        )
    }
    fn __gt__(&self, other: &Self) -> bool {
        self.0.partial_cmp(&other.0) == Some(Ordering::Greater)
    }
    fn __ge__(&self, other: &Self) -> bool {
        matches!(
            self.0.partial_cmp(&other.0),
            Some(Ordering::Greater | Ordering::Equal)
        )
    }

    // --- Number-specific ---

    fn max(&self, other: &Self) -> Self {
        Self(self.0.max(&other.0))
    }
    fn min(&self, other: &Self) -> Self {
        Self(self.0.min(&other.0))
    }
    fn clamp(&self, lo: &Self, hi: &Self) -> Self {
        Self(self.0.clamp(&lo.0, &hi.0))
    }

    // === Numeric elementary functions ===

    // --- Trigonometric ---
    fn sin(&self) -> Self {
        Self(self.0.sin())
    }
    fn cos(&self) -> Self {
        Self(self.0.cos())
    }
    fn tan(&self) -> Self {
        Self(self.0.tan())
    }
    fn cot(&self) -> Self {
        Self(self.0.cot())
    }
    fn sec(&self) -> Self {
        Self(self.0.sec())
    }
    fn csc(&self) -> Self {
        Self(self.0.csc())
    }

    // --- Inverse Trigonometric ---
    fn asin(&self) -> Self {
        Self(self.0.asin())
    }
    fn acos(&self) -> Self {
        Self(self.0.acos())
    }
    fn atan(&self) -> Self {
        Self(self.0.atan())
    }
    fn acot(&self) -> Self {
        Self(self.0.acot())
    }
    fn asec(&self) -> Self {
        Self(self.0.asec())
    }
    fn acsc(&self) -> Self {
        Self(self.0.acsc())
    }

    // --- Hyperbolic ---
    fn sinh(&self) -> Self {
        Self(self.0.sinh())
    }
    fn cosh(&self) -> Self {
        Self(self.0.cosh())
    }
    fn tanh(&self) -> Self {
        Self(self.0.tanh())
    }
    fn coth(&self) -> Self {
        Self(self.0.coth())
    }
    fn sech(&self) -> Self {
        Self(self.0.sech())
    }
    fn csch(&self) -> Self {
        Self(self.0.csch())
    }

    // --- Inverse Hyperbolic ---
    fn asinh(&self) -> Self {
        Self(self.0.asinh())
    }
    fn acosh(&self) -> Self {
        Self(self.0.acosh())
    }
    fn atanh(&self) -> Self {
        Self(self.0.atanh())
    }
    fn acoth(&self) -> Self {
        Self(self.0.acoth())
    }
    fn asech(&self) -> Self {
        Self(self.0.asech())
    }
    fn acsch(&self) -> Self {
        Self(self.0.acsch())
    }

    // --- Exponential & Logarithmic ---
    fn exp(&self) -> Self {
        Self(self.0.exp())
    }
    fn expm1(&self) -> Self {
        Self(self.0.expm1())
    }
    fn exp_neg(&self) -> Self {
        Self(self.0.exp_neg())
    }
    fn ln(&self) -> Self {
        Self(self.0.ln())
    }
    fn log1p(&self) -> Self {
        Self(self.0.log1p())
    }

    // --- Powers / Roots ---
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }
    fn cbrt(&self) -> Self {
        Self(self.0.cbrt())
    }

    // --- Basic Math ---
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn signum(&self) -> Self {
        Self(self.0.signum())
    }
    fn floor(&self) -> Self {
        Self(self.0.floor())
    }
    fn ceil(&self) -> Self {
        Self(self.0.ceil())
    }
    fn round(&self) -> Self {
        Self(self.0.round())
    }
    fn fract(&self) -> Self {
        Self(self.0.fract())
    }
    fn negate(&self) -> Self {
        Self(self.0.negate())
    }

    // --- Cardinal sine ---
    fn sinc(&self) -> Self {
        Self(self.0.sinc())
    }

    // --- Binary ---
    fn atan2(&self, x: &Self) -> Self {
        Self(self.0.atan2(&x.0))
    }
    fn log_base(&self, base: &Self) -> Self {
        Self(self.0.log_base(&base.0))
    }
    fn pow(&self, exp: &Self) -> Self {
        Self(self.0.pow(&exp.0))
    }

    // --- Properties ---
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
    fn is_neg_one(&self) -> bool {
        self.0.is_neg_one()
    }
    fn is_int(&self) -> bool {
        self.0.is_int()
    }
    fn is_int_ring(&self) -> bool {
        self.0.is_int_ring()
    }
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }
    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }
    fn to_float(&self) -> Self {
        Self(self.0.to_float())
    }
    fn approx_eq(&self, other: &Self, tol: &Self) -> bool {
        self.0.approx_eq_number(&other.0, &tol.0)
    }
    fn total_cmp(&self, other: &Self) -> i8 {
        match self.0.total_cmp(&other.0) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        }
    }
}

// ============================================================================
// PyCliffordNumber
// ============================================================================

#[cfg(feature = "clifford")]
#[pyclass(name = "CliffordNumber", skip_from_py_object)]
#[derive(Clone)]
struct PyCliffordNumber(CliffordNumber);

#[cfg(feature = "clifford")]
#[pymethods]
impl PyCliffordNumber {
    // --- Constructors ---

    #[staticmethod]
    fn zero(gens: &PyGeneratorSet) -> PyResult<Self> {
        CliffordNumber::zero(gens.0.clone())
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    #[staticmethod]
    fn scalar(gens: &PyGeneratorSet, value: &PyNumber) -> PyResult<Self> {
        CliffordNumber::scalar(gens.0.clone(), value.0.clone())
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    #[staticmethod]
    fn generator(gens: &PyGeneratorSet, index: u8) -> PyResult<Self> {
        CliffordNumber::generator(gens.0.clone(), index)
            .map(Self)
            .map_err(|e| map_error(&e))
    }

    // --- Display ---

    fn __repr__(&self) -> String {
        format!("CliffordNumber({})", self.0)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }

    // --- Accessors ---

    const fn blade_count(&self) -> usize {
        self.0.blade_count()
    }
    const fn n_generators(&self) -> usize {
        self.0.n_generators()
    }
    fn generator_set(&self) -> PyGeneratorSet {
        PyGeneratorSet(self.0.generator_set().clone())
    }
    fn coeff(&self, blade: usize) -> PyNumber {
        PyNumber(self.0.coeff(blade).clone().into())
    }
    fn set_coeff(&mut self, blade: usize, value: &PyNumber) -> PyResult<()> {
        self.0.set_coeff(
            blade,
            value.0.clone().into_real().map_err(|e| map_error(&e))?,
        );
        Ok(())
    }
    fn nonzero_blades(&self) -> Vec<(usize, PyNumber)> {
        self.0
            .nonzero_blades()
            .map(|(b, c)| (b, PyNumber(c.clone().into())))
            .collect()
    }
    // --- Geometric operations ---

    fn geometric_mul(&self, other: &Self) -> Self {
        Self(self.0.geometric_mul(&other.0))
    }
    fn outer_product(&self, other: &Self) -> Self {
        Self(self.0.outer_product(&other.0))
    }
    fn left_contraction(&self, other: &Self) -> Self {
        Self(self.0.left_contraction(&other.0))
    }
    fn scalar_product(&self, other: &Self) -> PyNumber {
        PyNumber(self.0.scalar_product(&other.0).into())
    }
    fn geometric_inverse(&self) -> Self {
        Self(self.0.geometric_inverse())
    }

    // --- Grade operations ---

    fn grade(&self, k: u32) -> Self {
        Self(self.0.grade(k))
    }
    fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }
    fn grade_involution(&self) -> Self {
        Self(self.0.grade_involution())
    }
    fn clifford_conjugate(&self) -> Self {
        Self(self.0.clifford_conjugate())
    }
    fn norm_sq(&self) -> PyNumber {
        PyNumber(self.0.norm_sq().into())
    }

    // --- Arithmetic ---

    fn __add__(&self, other: &Self) -> Self {
        Self(&self.0 + &other.0)
    }
    fn __sub__(&self, other: &Self) -> Self {
        Self(&self.0 - &other.0)
    }
    fn __mul__(&self, other: &Self) -> Self {
        Self(self.0.geometric_mul(&other.0))
    }
    fn __truediv__(&self, other: &Self) -> Self {
        Self(&self.0 / &other.0)
    }
    fn __neg__(&self) -> Self {
        Self(-&self.0)
    }
    fn __abs__(&self) -> Self {
        Self(self.0.abs())
    }
    fn __pow__(&self, exp: &Self, _mod: Option<&Self>) -> Self {
        Self(self.0.pow(&exp.0))
    }

    // --- Comparisons ---

    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __ne__(&self, other: &Self) -> bool {
        self.0 != other.0
    }

    // === Numeric elementary functions ===

    // --- Trigonometric ---
    fn sin(&self) -> Self {
        Self(self.0.sin())
    }
    fn cos(&self) -> Self {
        Self(self.0.cos())
    }
    fn tan(&self) -> Self {
        Self(self.0.tan())
    }
    fn cot(&self) -> Self {
        Self(self.0.cot())
    }
    fn sec(&self) -> Self {
        Self(self.0.sec())
    }
    fn csc(&self) -> Self {
        Self(self.0.csc())
    }

    // --- Inverse Trigonometric ---
    fn asin(&self) -> Self {
        Self(self.0.asin())
    }
    fn acos(&self) -> Self {
        Self(self.0.acos())
    }
    fn atan(&self) -> Self {
        Self(self.0.atan())
    }
    fn acot(&self) -> Self {
        Self(self.0.acot())
    }
    fn asec(&self) -> Self {
        Self(self.0.asec())
    }
    fn acsc(&self) -> Self {
        Self(self.0.acsc())
    }

    // --- Hyperbolic ---
    fn sinh(&self) -> Self {
        Self(self.0.sinh())
    }
    fn cosh(&self) -> Self {
        Self(self.0.cosh())
    }
    fn tanh(&self) -> Self {
        Self(self.0.tanh())
    }
    fn coth(&self) -> Self {
        Self(self.0.coth())
    }
    fn sech(&self) -> Self {
        Self(self.0.sech())
    }
    fn csch(&self) -> Self {
        Self(self.0.csch())
    }

    // --- Inverse Hyperbolic ---
    fn asinh(&self) -> Self {
        Self(self.0.asinh())
    }
    fn acosh(&self) -> Self {
        Self(self.0.acosh())
    }
    fn atanh(&self) -> Self {
        Self(self.0.atanh())
    }
    fn acoth(&self) -> Self {
        Self(self.0.acoth())
    }
    fn asech(&self) -> Self {
        Self(self.0.asech())
    }
    fn acsch(&self) -> Self {
        Self(self.0.acsch())
    }

    // --- Exponential & Logarithmic ---
    fn exp(&self) -> Self {
        Self(self.0.exp())
    }
    fn expm1(&self) -> Self {
        Self(self.0.expm1())
    }
    fn exp_neg(&self) -> Self {
        Self(self.0.exp_neg())
    }
    fn ln(&self) -> Self {
        Self(self.0.ln())
    }
    fn log1p(&self) -> Self {
        Self(self.0.log1p())
    }

    // --- Powers / Roots ---
    fn sqrt(&self) -> Self {
        Self(self.0.sqrt())
    }
    fn cbrt(&self) -> Self {
        Self(self.0.cbrt())
    }

    // --- Basic Math ---
    fn abs(&self) -> Self {
        Self(self.0.abs())
    }
    fn signum(&self) -> Self {
        Self(self.0.signum())
    }
    fn floor(&self) -> Self {
        Self(self.0.floor())
    }
    fn ceil(&self) -> Self {
        Self(self.0.ceil())
    }
    fn round(&self) -> Self {
        Self(self.0.round())
    }
    fn fract(&self) -> Self {
        Self(self.0.fract())
    }
    fn negate(&self) -> Self {
        Self(self.0.negate())
    }

    // --- Cardinal sine ---
    fn sinc(&self) -> Self {
        Self(self.0.sinc())
    }

    // --- Binary ---
    fn atan2(&self, x: &Self) -> Self {
        Self(self.0.atan2(&x.0))
    }
    fn log_base(&self, base: &Self) -> Self {
        Self(self.0.log_base(&base.0))
    }
    fn pow(&self, exp: &Self) -> Self {
        Self(self.0.pow(&exp.0))
    }

    // --- Properties ---
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_one(&self) -> bool {
        self.0.is_one()
    }
    fn is_neg_one(&self) -> bool {
        self.0.is_neg_one()
    }
    fn is_int(&self) -> bool {
        self.0.is_int()
    }
    fn is_int_ring(&self) -> bool {
        self.0.is_int_ring()
    }
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }
    fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
    fn is_positive(&self) -> bool {
        self.0.is_positive()
    }
    fn to_float(&self) -> Self {
        Self(self.0.to_float())
    }
    fn approx_eq(&self, other: &Self, tol: &Self) -> bool {
        self.0.approx_eq_number(&other.0, &tol.0)
    }
    fn total_cmp(&self, other: &Self) -> i8 {
        match self.0.total_cmp(&other.0) {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        }
    }
    fn max(&self, other: &Self) -> Self {
        Self(self.0.max(&other.0))
    }
    fn min(&self, other: &Self) -> Self {
        Self(self.0.min(&other.0))
    }
}

// ============================================================================
// PyGeneratorSet
// ============================================================================

#[cfg(feature = "clifford")]
#[pyclass(name = "GeneratorSet", skip_from_py_object)]
#[derive(Clone)]
struct PyGeneratorSet(GeneratorSet);

#[cfg(feature = "clifford")]
#[pymethods]
impl PyGeneratorSet {
    #[staticmethod]
    fn empty() -> Self {
        Self(GeneratorSet::empty())
    }

    #[staticmethod]
    fn cga() -> Self {
        Self(Cga::gens())
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.0)
    }
    fn __str__(&self) -> String {
        format!("{:?}", self.0)
    }
    const fn __len__(&self) -> usize {
        self.0.len()
    }
    fn __eq__(&self, other: &Self) -> bool {
        self.0 == other.0
    }
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.0.hash(&mut hasher);
        hasher.finish()
    }

    const fn len(&self) -> usize {
        self.0.len()
    }
    const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn id_at(&self, i: usize) -> u32 {
        self.0.id_at(i)
    }
    fn metric_at(&self, i: usize) -> i8 {
        self.0.metric_at(i)
    }
}

// ============================================================================
// Clifford free functions
// ============================================================================

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e1() -> PyCliffordNumber {
    PyCliffordNumber(Cga::e1())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e2() -> PyCliffordNumber {
    PyCliffordNumber(Cga::e2())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e3() -> PyCliffordNumber {
    PyCliffordNumber(Cga::e3())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e_plus() -> PyCliffordNumber {
    PyCliffordNumber(Cga::e_plus())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_e_minus() -> PyCliffordNumber {
    PyCliffordNumber(Cga::e_minus())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_orig() -> PyCliffordNumber {
    PyCliffordNumber(Cga::orig())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_inf() -> PyCliffordNumber {
    PyCliffordNumber(Cga::inf())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qi() -> PyCliffordNumber {
    PyCliffordNumber(Cga::qi())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qj() -> PyCliffordNumber {
    PyCliffordNumber(Cga::qj())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_qk() -> PyCliffordNumber {
    PyCliffordNumber(Cga::qk())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_pseudo3d() -> PyCliffordNumber {
    PyCliffordNumber(Cga::pseudo3d())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_pseudo5d() -> PyCliffordNumber {
    PyCliffordNumber(Cga::pseudo5d())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_eps() -> PyCliffordNumber {
    PyCliffordNumber(Cga::eps())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_ci() -> PyCliffordNumber {
    PyCliffordNumber(Cga::ci())
}

#[cfg(feature = "clifford")]
#[pyfunction]
fn py_sj() -> PyCliffordNumber {
    PyCliffordNumber(Cga::sj())
}

// ============================================================================
// Common free functions
// ============================================================================

#[pyfunction]
fn n(x: &Bound<'_, PyAny>) -> PyResult<PyNumber> {
    PyNumber::new(x)
}

#[pyfunction]
fn r(num: i64, den: i64) -> PyNumber {
    PyNumber(raw_r(num, den))
}

// ============================================================================
// Module registration
// ============================================================================

#[pymodule]
fn num_anafis_py(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyNumber>()?;
    m.add_function(wrap_pyfunction!(n, m)?)?;
    m.add_function(wrap_pyfunction!(r, m)?)?;

    #[cfg(feature = "clifford")]
    {
        m.add_class::<PyCliffordNumber>()?;
        m.add_class::<PyGeneratorSet>()?;
        m.add_function(wrap_pyfunction!(py_e1, m)?)?;
        m.add_function(wrap_pyfunction!(py_e2, m)?)?;
        m.add_function(wrap_pyfunction!(py_e3, m)?)?;
        m.add_function(wrap_pyfunction!(py_e_plus, m)?)?;
        m.add_function(wrap_pyfunction!(py_e_minus, m)?)?;
        m.add_function(wrap_pyfunction!(py_orig, m)?)?;
        m.add_function(wrap_pyfunction!(py_inf, m)?)?;
        m.add_function(wrap_pyfunction!(py_qi, m)?)?;
        m.add_function(wrap_pyfunction!(py_qj, m)?)?;
        m.add_function(wrap_pyfunction!(py_qk, m)?)?;
        m.add_function(wrap_pyfunction!(py_pseudo3d, m)?)?;
        m.add_function(wrap_pyfunction!(py_pseudo5d, m)?)?;
        m.add_function(wrap_pyfunction!(py_eps, m)?)?;
        m.add_function(wrap_pyfunction!(py_ci, m)?)?;
        m.add_function(wrap_pyfunction!(py_sj, m)?)?;
    }

    Ok(())
}

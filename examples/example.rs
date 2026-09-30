//! Comprehensive example showcasing *everything* exported by `num_anafis`!
//!
//! Run with: cargo run --example example --features "clifford"

#![expect(
    clippy::print_stdout,
    clippy::use_debug,
    clippy::uninlined_format_args,
    clippy::non_ascii_literal,
    clippy::too_many_lines,
    reason = "Examples must print and debug-print values; convenience over strict linting"
)]

use num_anafis::{F64Ext, FloatType, IntType, NumAnafisError, Numeric, RationalType, c, i, n, r};

#[cfg(feature = "clifford")]
use num_anafis::{Cga, CliffordNumber, FastClifford, GeneratorSet};

fn title(t: &str) {
    println!("\n━━━ {} ━━━", t);
}

#[cfg(feature = "clifford")]
fn put(label: &str, val: &CliffordNumber) {
    println!("  {:<25} = {}", label, val);
}

fn main() {
    println!("============================================================");
    println!("  Welcome to num-anafis: Numeric CAS Foundation            ");
    println!("============================================================");

    // ==================================================================
    // 1. Primitive Types & Extension Traits
    // ==================================================================
    title("1. Types & Traits (IntType, FloatType, RationalType, F64Ext)");
    let my_float: FloatType = core::f64::consts::PI;

    println!(
        "  FloatType ({})         = {}",
        core::any::type_name::<FloatType>(),
        my_float
    );
    println!(
        "  my_float.lgamma()       = {} (via F64Ext)",
        my_float.lgamma()
    );

    let my_int: IntType = IntType::from(42);

    println!(
        "  IntType ({})           = {}",
        core::any::type_name::<IntType>(),
        my_int
    );

    println!(
        "  RationalType ({})      = fractional backend",
        core::any::type_name::<RationalType>()
    );

    // ==================================================================
    // 2. Numbers & Conversions (4 Tiers: Int -> Rational -> Float -> Complex)
    // ==================================================================
    title("2. 4-Tier Numbers (Int, Rational, Float, Complex)");
    let num1 = n(42);
    let num2 = r(22, 7);
    let num3: num_anafis::Number = 100.into();
    let num4: num_anafis::Number = core::f64::consts::PI.into();
    let num_c = c(3, 4);
    let num_f = n(2.5);
    let num_im = 5 * i();
    let num_i = i();

    println!("  Integer 'n(42)'       = {}", num1);
    println!("  Rational 'r(22, 7)'   = {}", num2);
    println!("  Float 'n(2.5)'        = {}", num_f);
    println!("  100.into()     = {}", num3);
    println!("  3.14.into()    = {}", num4);
    println!("  Complex 'c(3, 4)'     = {}", num_c);
    println!("  Imaginary '5 * i()'     = {}", num_im);
    println!("  Imaginary unit 'i()'  = {}", num_i);
    println!("  sqrt(-1) in CAS       = {}", n(-1).sqrt());
    println!("  ln(-1) in CAS         = {}", n(-1).ln());
    println!("  Numeric trait (exp)   = {}", num1.exp());

    // ==================================================================
    // 3. Error Handling
    // ==================================================================
    title("3. Error Handling (NumAnafisError)");
    #[cfg(feature = "clifford")]
    {
        let zero_mv = CliffordNumber::scalar(GeneratorSet::empty(), n(0))
            .expect("empty algebra admits real scalars");
        println!(
            "  Inverse of 0          = {:?}",
            zero_mv.geometric_inverse()
        );
    }
    #[cfg(not(feature = "clifford"))]
    println!("  (clifford feature disabled — skipping geometric algebra examples)");
    let err_example = NumAnafisError::InvalidGeneratorMetric { id: 256, metric: 2 };
    println!("  Demonstrating Error   = {:?}", err_example);

    // ==================================================================
    // Clifford-only sections (feature gate)
    // ==================================================================
    #[cfg(feature = "clifford")]
    {
        // ==============================================================
        // 4. Generator Sets (GeneratorSet, cga_gens)
        // ==============================================================
        title("4. Generator Sets");
        let empty_gen = GeneratorSet::empty();
        let cga = Cga::gens();
        println!("  Empty GeneratorSet    = {:?}", empty_gen);
        println!("  CGA GeneratorSet      = {} dimensions active", cga.len());

        // ==============================================================
        // 5. Convenience Constructors (e1..e3, ci, sj, eps, etc.)
        // ==============================================================
        title("5. Built-in Geometric Algebra Generators");
        put("Complex 'ci' (ci² = -1)", &Cga::ci());
        put("Split-Complex 'sj' (sj² = +1)", &Cga::sj());
        put("Dual 'ε' (ε² = 0)", &Cga::eps());
        put("Quaternions 'qi'", &Cga::qi());
        put("Quaternions 'qj'", &Cga::qj());
        put("Quaternions 'qk'", &Cga::qk());

        put("e1 (e1² = +1)", &Cga::e1());
        put("e2 (e2² = +1)", &Cga::e2());
        put("e3 (e3² = +1)", &Cga::e3());

        put("3D Pseudoscalar", &Cga::pseudo3d());
        put("5D Pseudoscalar", &Cga::pseudo5d());

        // ==============================================================
        // 6. Conformal Geometric Algebra (CGA) Identifiers
        // ==============================================================
        title("6. Conformal Geometric Algebra (CGA) - Cl(4,1)");
        put("Extra positive 'e_plus'", &Cga::e_plus());
        put("Extra negative 'e_minus'", &Cga::e_minus());
        put("Origin 'orig'", &Cga::orig());
        put("Infinity 'inf'", &Cga::inf());

        put("orig * inf", &(&Cga::orig() * &Cga::inf()));

        // ==============================================================
        // 7. Dynamic Multivectors (CliffordNumber)
        // ==============================================================
        title("7. Dynamic Multivectors & Spectral Math (CliffordNumber)");
        let mut mv = &Cga::e1() * 3_i64 + &(&Cga::e1() * &Cga::e2()) * 4_i64;
        mv = &mv + &Cga::pseudo3d();
        put("Multivector 'M'", &mv);

        put("sin(M)", &mv.sin());
        put("exp(M)", &mv.exp());

        // ==============================================================
        // 8. FastClifford
        // ==============================================================
        title("8. FastClifford: Zero-Allocation Compile-Time Generic Math");
        let mut fc_a = FastClifford::<3, 0, 0>::zero();
        *fc_a
            .coeffs_mut_slice()
            .get_mut(0)
            .expect("Cl(3,0) includes a scalar blade") = num_anafis::Real::from_int(2);
        *fc_a
            .coeffs_mut_slice()
            .get_mut(1)
            .expect("Cl(3,0) includes its first generator") = num_anafis::Real::from_int(3);

        let mut fc_b = FastClifford::<3, 0, 0>::zero();
        *fc_b
            .coeffs_mut_slice()
            .get_mut(1)
            .expect("Cl(3,0) includes its first generator") = num_anafis::Real::from_int(5);

        println!("  FastClifford A        = {}", fc_a);
        println!("  FastClifford B        = {}", fc_b);
        println!("  A + B                 = {}", fc_a.clone() + fc_b.clone());
        println!("  A * B                 = {}", fc_a.clone() * fc_b);

        println!("  exp(A)                = {}", fc_a.exp());
        println!("  sin(A)                = {}", fc_a.sin());
    }

    println!("\n num-anafis demonstration complete!");
}

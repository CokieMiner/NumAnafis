//! Binary batch evaluator for the precision CLI's primitive extension traits.

use core::error::Error;

use std::{
    env,
    io::{self, BufReader, BufWriter, Read, Write},
};

use num_anafis::{F32Ext, F64Ext};

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = env::args().skip(1);
    let width = arguments.next().ok_or("missing width: 32 or 64")?;
    let function = arguments.next().ok_or("missing function")?;
    if arguments.next().is_some() || !matches!(width.as_str(), "32" | "64") {
        return Err("expected width and function".into());
    }
    let mut input = BufReader::new(io::stdin().lock());
    let mut output = BufWriter::new(io::stdout().lock());
    let mut bytes = [0_u8; 32];
    loop {
        if input.read(&mut bytes[..1])? == 0 {
            break;
        }
        input.read_exact(&mut bytes[1..])?;
        let mut row = [0_u64; 4];
        for (word, chunk) in row.iter_mut().zip(bytes.as_chunks::<8>().0) {
            *word = u64::from_le_bytes(*chunk);
        }
        let value = if width == "32" {
            evaluate32(&function, row)?
        } else {
            evaluate64(&function, row)?
        };
        output.write_all(&value.to_le_bytes())?;
        output.flush()?;
    }
    drop(input);
    output.flush()?;
    Ok(())
}

fn evaluate32(function: &str, row: [u64; 4]) -> Result<u64, Box<dyn Error>> {
    let x = f32::from_bits(u32::try_from(row[0])?);
    let result = match function {
        "sin" => F32Ext::sin(x),
        "cos" => F32Ext::cos(x),
        "tan" => F32Ext::tan(x),
        "asin" => F32Ext::asin(x),
        "acos" => F32Ext::acos(x),
        "atan" => F32Ext::atan(x),
        "sinh" => F32Ext::sinh(x),
        "cosh" => F32Ext::cosh(x),
        "tanh" => F32Ext::tanh(x),
        "asinh" => F32Ext::asinh(x),
        "acosh" => F32Ext::acosh(x),
        "atanh" => F32Ext::atanh(x),
        "exp" => F32Ext::exp(x),
        "exp_m1" => F32Ext::exp_m1(x),
        "ln" => F32Ext::ln(x),
        "ln_1p" => F32Ext::ln_1p(x),
        "sqrt" => F32Ext::sqrt(x),
        "cbrt" => F32Ext::cbrt(x),
        "floor" => F32Ext::floor(x),
        "ceil" => F32Ext::ceil(x),
        "round" => F32Ext::round(x),
        "trunc" => F32Ext::trunc(x),
        "fract" => F32Ext::fract(x),
        "gamma" => F32Ext::gamma(x),
        "lgamma" => F32Ext::lgamma(x),
        "digamma" => F32Ext::digamma(x),
        "trigamma" => F32Ext::trigamma(x),
        "tetragamma" => F32Ext::tetragamma(x),
        "erf" => F32Ext::erf(x),
        "erfc" => F32Ext::erfc(x),
        "zeta" => F32Ext::zeta(x),
        "elliptic_k" => F32Ext::elliptic_k(x),
        "elliptic_e" => F32Ext::elliptic_e(x),
        "lambert_w0" => F32Ext::lambert_w0(x),
        "lambert_wm1" => F32Ext::lambert_wm1(x),
        "atan2" => F32Ext::atan2(x, f32::from_bits(u32::try_from(row[1])?)),
        "powf" => F32Ext::powf(x, f32::from_bits(u32::try_from(row[1])?)),
        "beta" => F32Ext::beta(x, f32::from_bits(u32::try_from(row[1])?)),
        "polygamma" => {
            F32Ext::polygamma(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?)
        }
        "zeta_deriv" => {
            F32Ext::zeta_deriv(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?)
        }
        "bessel_j" => F32Ext::bessel_j(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?),
        "bessel_y" => F32Ext::bessel_y(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?),
        "bessel_i" => F32Ext::bessel_i(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?),
        "bessel_k" => F32Ext::bessel_k(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?),
        "hermite" => F32Ext::hermite(x, i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?),
        "assoc_legendre" => F32Ext::assoc_legendre(
            x,
            i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?,
            i32::try_from(i64::from_le_bytes(row[2].to_le_bytes()))?,
        ),
        "spherical_harmonic" => F32Ext::spherical_harmonic(
            x,
            i32::try_from(i64::from_le_bytes(row[1].to_le_bytes()))?,
            i32::try_from(i64::from_le_bytes(row[2].to_le_bytes()))?,
            f32::from_bits(u32::try_from(row[3])?),
        ),
        _ => return Err(format!("unknown function: {function}").into()),
    };
    Ok(u64::from(result.to_bits()))
}

fn evaluate64(function: &str, row: [u64; 4]) -> Result<u64, Box<dyn Error>> {
    let x = f64::from_bits(row[0]);
    let result = match function {
        "sin" => F64Ext::sin(x),
        "cos" => F64Ext::cos(x),
        "tan" => F64Ext::tan(x),
        "asin" => F64Ext::asin(x),
        "acos" => F64Ext::acos(x),
        "atan" => F64Ext::atan(x),
        "sinh" => F64Ext::sinh(x),
        "cosh" => F64Ext::cosh(x),
        "tanh" => F64Ext::tanh(x),
        "asinh" => F64Ext::asinh(x),
        "acosh" => F64Ext::acosh(x),
        "atanh" => F64Ext::atanh(x),
        "exp" => F64Ext::exp(x),
        "exp_m1" => F64Ext::exp_m1(x),
        "ln" => F64Ext::ln(x),
        "ln_1p" => F64Ext::ln_1p(x),
        "sqrt" => F64Ext::sqrt(x),
        "cbrt" => F64Ext::cbrt(x),
        "floor" => F64Ext::floor(x),
        "ceil" => F64Ext::ceil(x),
        "round" => F64Ext::round(x),
        "trunc" => F64Ext::trunc(x),
        "fract" => F64Ext::fract(x),
        "gamma" => F64Ext::gamma(x),
        "lgamma" => F64Ext::lgamma(x),
        "digamma" => F64Ext::digamma(x),
        "trigamma" => F64Ext::trigamma(x),
        "tetragamma" => F64Ext::tetragamma(x),
        "erf" => F64Ext::erf(x),
        "erfc" => F64Ext::erfc(x),
        "zeta" => F64Ext::zeta(x),
        "elliptic_k" => F64Ext::elliptic_k(x),
        "elliptic_e" => F64Ext::elliptic_e(x),
        "lambert_w0" => F64Ext::lambert_w0(x),
        "lambert_wm1" => F64Ext::lambert_wm1(x),
        "atan2" => F64Ext::atan2(x, f64::from_bits(row[1])),
        "powf" => F64Ext::powf(x, f64::from_bits(row[1])),
        "beta" => F64Ext::beta(x, f64::from_bits(row[1])),
        "polygamma" => F64Ext::polygamma(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "zeta_deriv" => F64Ext::zeta_deriv(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "bessel_j" => F64Ext::bessel_j(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "bessel_y" => F64Ext::bessel_y(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "bessel_i" => F64Ext::bessel_i(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "bessel_k" => F64Ext::bessel_k(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "hermite" => F64Ext::hermite(x, i64::from_le_bytes(row[1].to_le_bytes())),
        "assoc_legendre" => F64Ext::assoc_legendre(
            x,
            i64::from_le_bytes(row[1].to_le_bytes()),
            i64::from_le_bytes(row[2].to_le_bytes()),
        ),
        "spherical_harmonic" => F64Ext::spherical_harmonic(
            x,
            i64::from_le_bytes(row[1].to_le_bytes()),
            i64::from_le_bytes(row[2].to_le_bytes()),
            f64::from_bits(row[3]),
        ),
        _ => return Err(format!("unknown function: {function}").into()),
    };
    Ok(result.to_bits())
}

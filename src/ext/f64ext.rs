//! Extension trait and mathematical functions for `f64`.

#[cfg(not(feature = "std"))]
use libm::{
    acos, acosh, asin, asinh, atan, atan2, atanh, cbrt, ceil, cos, cosh, exp, expm1, floor, log,
    log1p, pow, round, sin, sinh, sqrt, tan, tanh, trunc,
};

use crate::special::SpecialFunc;

/// Mathematical extensions and special functions for [`f64`].
pub trait F64Ext {
    /// Sine.
    fn sin(self) -> f64;
    /// Cosine.
    fn cos(self) -> f64;
    /// Tangent.
    fn tan(self) -> f64;
    /// Inverse sine.
    fn asin(self) -> f64;
    /// Inverse cosine.
    fn acos(self) -> f64;
    /// Inverse tangent.
    fn atan(self) -> f64;
    /// Four-quadrant inverse tangent.
    fn atan2(self, other: f64) -> f64;

    /// Hyperbolic sine.
    fn sinh(self) -> f64;
    /// Hyperbolic cosine.
    fn cosh(self) -> f64;
    /// Hyperbolic tangent.
    fn tanh(self) -> f64;
    /// Inverse hyperbolic sine.
    fn asinh(self) -> f64;
    /// Inverse hyperbolic cosine.
    fn acosh(self) -> f64;
    /// Inverse hyperbolic tangent.
    fn atanh(self) -> f64;

    /// Exponential $e^x$.
    fn exp(self) -> f64;
    /// Exponential minus one $e^x - 1$.
    fn exp_m1(self) -> f64;
    /// Natural logarithm $\ln(x)$.
    fn ln(self) -> f64;
    /// Natural logarithm of $1 + x$.
    fn ln_1p(self) -> f64;
    /// Square root $\sqrt{x}$.
    fn sqrt(self) -> f64;
    /// Cube root $\sqrt\[3\]{x}$.
    fn cbrt(self) -> f64;
    /// Float power $x^y$.
    fn powf(self, n: f64) -> f64;
    /// Largest integer less than or equal to `self`.
    fn floor(self) -> f64;
    /// Smallest integer greater than or equal to `self`.
    fn ceil(self) -> f64;
    /// Nearest integer, rounding half away from zero.
    fn round(self) -> f64;
    /// Integer part, truncating toward zero.
    fn trunc(self) -> f64;
    /// Fractional part (`self` minus its truncation).
    fn fract(self) -> f64;

    /// Gamma function $\Gamma(x)$.
    fn gamma(self) -> f64;
    /// Log-gamma function $\ln|\Gamma(x)|$.
    fn lgamma(self) -> f64;
    /// Digamma function $\psi(x) = \Gamma'(x)/\Gamma(x)$.
    fn digamma(self) -> f64;
    /// Trigamma function $\psi^{(1)}(x)$.
    fn trigamma(self) -> f64;
    /// Tetragamma function $\psi^{(2)}(x)$.
    fn tetragamma(self) -> f64;
    /// Polygamma function $\psi^{(n)}(x)$.
    fn polygamma(self, n: i64) -> f64;
    /// Error function $\operatorname{erf}(x)$.
    fn erf(self) -> f64;
    /// Complementary error function $\operatorname{erfc}(x)$.
    fn erfc(self) -> f64;
    /// Riemann zeta function $\zeta(x)$.
    fn zeta(self) -> f64;
    /// $n$-th derivative of the Riemann zeta function $\zeta^{(n)}(x)$.
    fn zeta_deriv(self, n: i64) -> f64;
    /// Bessel function of the first kind $`J_n(x)`$.
    fn bessel_j(self, n: i64) -> f64;
    /// Bessel function of the second kind $`Y_n(x)`$.
    fn bessel_y(self, n: i64) -> f64;
    /// Modified Bessel function of the first kind $`I_n(x)`$.
    fn bessel_i(self, n: i64) -> f64;
    /// Modified Bessel function of the second kind $`K_n(x)`$.
    fn bessel_k(self, n: i64) -> f64;
    /// Beta function $\mathrm{B}(x, y)$.
    fn beta(self, b: f64) -> f64;
    /// Complete elliptic integral of the first kind $K(k)$.
    fn elliptic_k(self) -> f64;
    /// Complete elliptic integral of the second kind $E(k)$.
    fn elliptic_e(self) -> f64;
    /// Physicists' Hermite polynomial $`H_n(x)`$.
    fn hermite(self, n: i64) -> f64;
    /// Principal branch of the Lambert W function $`W_0(x)`$.
    fn lambert_w0(self) -> f64;
    /// Secondary real branch of the Lambert W function $W_{-1}(x)$.
    fn lambert_wm1(self) -> f64;
    /// Associated Legendre polynomial $`P_l^m(x)`$.
    fn assoc_legendre(self, l: i64, m: i64) -> f64;
    /// Spherical harmonic $`Y_l^m(\theta`, \phi)$.
    fn spherical_harmonic(self, l: i64, m: i64, phi: f64) -> f64;
}

impl F64Ext for f64 {
    #[inline]
    fn sin(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.sin()
        }
        #[cfg(not(feature = "std"))]
        {
            sin(self)
        }
    }

    #[inline]
    fn cos(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.cos()
        }
        #[cfg(not(feature = "std"))]
        {
            cos(self)
        }
    }

    #[inline]
    fn tan(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.tan()
        }
        #[cfg(not(feature = "std"))]
        {
            tan(self)
        }
    }

    #[inline]
    fn asin(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.asin()
        }
        #[cfg(not(feature = "std"))]
        {
            asin(self)
        }
    }

    #[inline]
    fn acos(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.acos()
        }
        #[cfg(not(feature = "std"))]
        {
            acos(self)
        }
    }

    #[inline]
    fn atan(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.atan()
        }
        #[cfg(not(feature = "std"))]
        {
            atan(self)
        }
    }

    #[inline]
    fn atan2(self, other: f64) -> f64 {
        #[cfg(feature = "std")]
        {
            self.atan2(other)
        }
        #[cfg(not(feature = "std"))]
        {
            atan2(self, other)
        }
    }

    #[inline]
    fn sinh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.sinh()
        }
        #[cfg(not(feature = "std"))]
        {
            sinh(self)
        }
    }

    #[inline]
    fn cosh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.cosh()
        }
        #[cfg(not(feature = "std"))]
        {
            cosh(self)
        }
    }

    #[inline]
    fn tanh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.tanh()
        }
        #[cfg(not(feature = "std"))]
        {
            tanh(self)
        }
    }

    #[inline]
    fn asinh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.asinh()
        }
        #[cfg(not(feature = "std"))]
        {
            asinh(self)
        }
    }

    #[inline]
    fn acosh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.acosh()
        }
        #[cfg(not(feature = "std"))]
        {
            acosh(self)
        }
    }

    #[inline]
    fn atanh(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.atanh()
        }
        #[cfg(not(feature = "std"))]
        {
            atanh(self)
        }
    }

    #[inline]
    fn exp(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.exp()
        }
        #[cfg(not(feature = "std"))]
        {
            exp(self)
        }
    }

    #[inline]
    fn exp_m1(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.exp_m1()
        }
        #[cfg(not(feature = "std"))]
        {
            expm1(self)
        }
    }

    #[inline]
    fn ln(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.ln()
        }
        #[cfg(not(feature = "std"))]
        {
            log(self)
        }
    }

    #[inline]
    fn ln_1p(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.ln_1p()
        }
        #[cfg(not(feature = "std"))]
        {
            log1p(self)
        }
    }

    #[inline]
    fn sqrt(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.sqrt()
        }
        #[cfg(not(feature = "std"))]
        {
            sqrt(self)
        }
    }

    #[inline]
    fn cbrt(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.cbrt()
        }
        #[cfg(not(feature = "std"))]
        {
            cbrt(self)
        }
    }

    #[inline]
    fn powf(self, n: f64) -> f64 {
        #[cfg(feature = "std")]
        {
            self.powf(n)
        }
        #[cfg(not(feature = "std"))]
        {
            pow(self, n)
        }
    }

    #[inline]
    fn floor(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.floor()
        }
        #[cfg(not(feature = "std"))]
        {
            floor(self)
        }
    }

    #[inline]
    fn ceil(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.ceil()
        }
        #[cfg(not(feature = "std"))]
        {
            ceil(self)
        }
    }

    #[inline]
    fn round(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.round()
        }
        #[cfg(not(feature = "std"))]
        {
            round(self)
        }
    }

    #[inline]
    fn trunc(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.trunc()
        }
        #[cfg(not(feature = "std"))]
        {
            trunc(self)
        }
    }

    #[inline]
    fn fract(self) -> f64 {
        #[cfg(feature = "std")]
        {
            self.fract()
        }
        // libm provides no fract: inherent is defined as self - self.trunc(),
        // reproduced exactly (including NaN, infinite, and signed-zero edges).
        #[cfg(not(feature = "std"))]
        {
            self - trunc(self)
        }
    }

    #[inline]
    fn gamma(self) -> f64 {
        SpecialFunc::gamma(self)
    }

    #[inline]
    fn lgamma(self) -> f64 {
        SpecialFunc::lgamma(self)
    }

    #[inline]
    fn digamma(self) -> f64 {
        SpecialFunc::digamma(self)
    }

    #[inline]
    fn trigamma(self) -> f64 {
        SpecialFunc::trigamma(self)
    }

    #[inline]
    fn tetragamma(self) -> f64 {
        SpecialFunc::tetragamma(self)
    }

    #[inline]
    fn polygamma(self, n: i64) -> f64 {
        SpecialFunc::polygamma_n(n, self)
    }

    #[inline]
    fn erf(self) -> f64 {
        SpecialFunc::erf(self)
    }

    #[inline]
    fn erfc(self) -> f64 {
        SpecialFunc::erfc(self)
    }

    #[inline]
    fn zeta(self) -> f64 {
        SpecialFunc::zeta(self)
    }

    #[inline]
    fn zeta_deriv(self, n: i64) -> f64 {
        SpecialFunc::zeta_deriv(n, self)
    }

    #[inline]
    fn bessel_j(self, n: i64) -> f64 {
        SpecialFunc::besselj(n, self)
    }

    #[inline]
    fn bessel_y(self, n: i64) -> f64 {
        SpecialFunc::bessely(n, self)
    }

    #[inline]
    fn bessel_i(self, n: i64) -> f64 {
        SpecialFunc::besseli(n, self)
    }

    #[inline]
    fn bessel_k(self, n: i64) -> f64 {
        SpecialFunc::besselk(n, self)
    }

    #[inline]
    fn beta(self, b: f64) -> f64 {
        SpecialFunc::beta(self, b)
    }

    #[inline]
    fn elliptic_k(self) -> f64 {
        SpecialFunc::elliptic_k(self)
    }

    #[inline]
    fn elliptic_e(self) -> f64 {
        SpecialFunc::elliptic_e(self)
    }

    #[inline]
    fn hermite(self, n: i64) -> f64 {
        SpecialFunc::hermite(n, self)
    }

    #[inline]
    fn lambert_w0(self) -> f64 {
        SpecialFunc::lambertw0(self)
    }

    #[inline]
    fn lambert_wm1(self) -> f64 {
        SpecialFunc::lambertwm1(self)
    }

    #[inline]
    fn assoc_legendre(self, l: i64, m: i64) -> f64 {
        SpecialFunc::assoc_legendre(l, m, self)
    }

    #[inline]
    fn spherical_harmonic(self, l: i64, m: i64, phi: f64) -> f64 {
        SpecialFunc::spherical_harmonic(l, m, self, phi)
    }
}

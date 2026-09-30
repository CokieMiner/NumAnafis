//! Extension trait and mathematical functions for `f32`.

#[cfg(not(feature = "std"))]
use libm::{
    acosf, acoshf, asinf, asinhf, atan2f, atanf, atanhf, cbrtf, ceilf, cosf, coshf, expf, expm1f,
    floorf, log1pf, logf, powf, roundf, sinf, sinhf, sqrtf, tanf, tanhf, truncf,
};

use crate::special::SpecialFunc;

/// Mathematical extensions and special functions for [`f32`].
pub trait F32Ext {
    /// Sine.
    fn sin(self) -> f32;
    /// Cosine.
    fn cos(self) -> f32;
    /// Tangent.
    fn tan(self) -> f32;
    /// Inverse sine.
    fn asin(self) -> f32;
    /// Inverse cosine.
    fn acos(self) -> f32;
    /// Inverse tangent.
    fn atan(self) -> f32;
    /// Four-quadrant inverse tangent.
    fn atan2(self, other: f32) -> f32;

    /// Hyperbolic sine.
    fn sinh(self) -> f32;
    /// Hyperbolic cosine.
    fn cosh(self) -> f32;
    /// Hyperbolic tangent.
    fn tanh(self) -> f32;
    /// Inverse hyperbolic sine.
    fn asinh(self) -> f32;
    /// Inverse hyperbolic cosine.
    fn acosh(self) -> f32;
    /// Inverse hyperbolic tangent.
    fn atanh(self) -> f32;

    /// Exponential $e^x$.
    fn exp(self) -> f32;
    /// Exponential minus one $e^x - 1$.
    fn exp_m1(self) -> f32;
    /// Natural logarithm $\ln(x)$.
    fn ln(self) -> f32;
    /// Natural logarithm of $1 + x$.
    fn ln_1p(self) -> f32;
    /// Square root $\sqrt{x}$.
    fn sqrt(self) -> f32;
    /// Cube root $\sqrt\[3\]{x}$.
    fn cbrt(self) -> f32;
    /// Float power $x^y$.
    fn powf(self, n: f32) -> f32;
    /// Largest integer less than or equal to `self`.
    fn floor(self) -> f32;
    /// Smallest integer greater than or equal to `self`.
    fn ceil(self) -> f32;
    /// Nearest integer, rounding half away from zero.
    fn round(self) -> f32;
    /// Integer part, truncating toward zero.
    fn trunc(self) -> f32;
    /// Fractional part (`self` minus its truncation).
    fn fract(self) -> f32;

    /// Gamma function $\Gamma(x)$.
    fn gamma(self) -> f32;
    /// Log-gamma function $\ln|\Gamma(x)|$.
    fn lgamma(self) -> f32;
    /// Digamma function $\psi(x) = \Gamma'(x)/\Gamma(x)$.
    fn digamma(self) -> f32;
    /// Trigamma function $\psi^{(1)}(x)$.
    fn trigamma(self) -> f32;
    /// Tetragamma function $\psi^{(2)}(x)$.
    fn tetragamma(self) -> f32;
    /// Polygamma function $\psi^{(n)}(x)$.
    fn polygamma(self, n: i32) -> f32;
    /// Error function $\operatorname{erf}(x)$.
    fn erf(self) -> f32;
    /// Complementary error function $\operatorname{erfc}(x)$.
    fn erfc(self) -> f32;
    /// Riemann zeta function $\zeta(x)$.
    fn zeta(self) -> f32;
    /// $n$-th derivative of the Riemann zeta function $\zeta^{(n)}(x)$.
    fn zeta_deriv(self, n: i32) -> f32;
    /// Bessel function of the first kind $`J_n(x)`$.
    fn bessel_j(self, n: i32) -> f32;
    /// Bessel function of the second kind $`Y_n(x)`$.
    fn bessel_y(self, n: i32) -> f32;
    /// Modified Bessel function of the first kind $`I_n(x)`$.
    fn bessel_i(self, n: i32) -> f32;
    /// Modified Bessel function of the second kind $`K_n(x)`$.
    fn bessel_k(self, n: i32) -> f32;
    /// Beta function $\mathrm{B}(x, y)$.
    fn beta(self, b: f32) -> f32;
    /// Complete elliptic integral of the first kind $K(k)$.
    fn elliptic_k(self) -> f32;
    /// Complete elliptic integral of the second kind $E(k)$.
    fn elliptic_e(self) -> f32;
    /// Physicists' Hermite polynomial $`H_n(x)`$.
    fn hermite(self, n: i32) -> f32;
    /// Principal branch of the Lambert W function $`W_0(x)`$.
    fn lambert_w0(self) -> f32;
    /// Secondary real branch of the Lambert W function $W_{-1}(x)$.
    fn lambert_wm1(self) -> f32;
    /// Associated Legendre polynomial $`P_l^m(x)`$.
    fn assoc_legendre(self, l: i32, m: i32) -> f32;
    /// Spherical harmonic $`Y_l^m(\theta`, \phi)$.
    fn spherical_harmonic(self, l: i32, m: i32, phi: f32) -> f32;
}

impl F32Ext for f32 {
    #[inline]
    fn sin(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.sin()
        }
        #[cfg(not(feature = "std"))]
        {
            sinf(self)
        }
    }

    #[inline]
    fn cos(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.cos()
        }
        #[cfg(not(feature = "std"))]
        {
            cosf(self)
        }
    }

    #[inline]
    fn tan(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.tan()
        }
        #[cfg(not(feature = "std"))]
        {
            tanf(self)
        }
    }

    #[inline]
    fn asin(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.asin()
        }
        #[cfg(not(feature = "std"))]
        {
            asinf(self)
        }
    }

    #[inline]
    fn acos(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.acos()
        }
        #[cfg(not(feature = "std"))]
        {
            acosf(self)
        }
    }

    #[inline]
    fn atan(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.atan()
        }
        #[cfg(not(feature = "std"))]
        {
            atanf(self)
        }
    }

    #[inline]
    fn atan2(self, other: f32) -> f32 {
        #[cfg(feature = "std")]
        {
            self.atan2(other)
        }
        #[cfg(not(feature = "std"))]
        {
            atan2f(self, other)
        }
    }

    #[inline]
    fn sinh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.sinh()
        }
        #[cfg(not(feature = "std"))]
        {
            sinhf(self)
        }
    }

    #[inline]
    fn cosh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.cosh()
        }
        #[cfg(not(feature = "std"))]
        {
            coshf(self)
        }
    }

    #[inline]
    fn tanh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.tanh()
        }
        #[cfg(not(feature = "std"))]
        {
            tanhf(self)
        }
    }

    #[inline]
    fn asinh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.asinh()
        }
        #[cfg(not(feature = "std"))]
        {
            asinhf(self)
        }
    }

    #[inline]
    fn acosh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.acosh()
        }
        #[cfg(not(feature = "std"))]
        {
            acoshf(self)
        }
    }

    #[inline]
    fn atanh(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.atanh()
        }
        #[cfg(not(feature = "std"))]
        {
            atanhf(self)
        }
    }

    #[inline]
    fn exp(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.exp()
        }
        #[cfg(not(feature = "std"))]
        {
            expf(self)
        }
    }

    #[inline]
    fn exp_m1(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.exp_m1()
        }
        #[cfg(not(feature = "std"))]
        {
            expm1f(self)
        }
    }

    #[inline]
    fn ln(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.ln()
        }
        #[cfg(not(feature = "std"))]
        {
            logf(self)
        }
    }

    #[inline]
    fn ln_1p(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.ln_1p()
        }
        #[cfg(not(feature = "std"))]
        {
            log1pf(self)
        }
    }

    #[inline]
    fn sqrt(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.sqrt()
        }
        #[cfg(not(feature = "std"))]
        {
            sqrtf(self)
        }
    }

    #[inline]
    fn cbrt(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.cbrt()
        }
        #[cfg(not(feature = "std"))]
        {
            cbrtf(self)
        }
    }

    #[inline]
    fn powf(self, n: f32) -> f32 {
        #[cfg(feature = "std")]
        {
            self.powf(n)
        }
        #[cfg(not(feature = "std"))]
        {
            powf(self, n)
        }
    }

    #[inline]
    fn floor(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.floor()
        }
        #[cfg(not(feature = "std"))]
        {
            floorf(self)
        }
    }

    #[inline]
    fn ceil(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.ceil()
        }
        #[cfg(not(feature = "std"))]
        {
            ceilf(self)
        }
    }

    #[inline]
    fn round(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.round()
        }
        #[cfg(not(feature = "std"))]
        {
            roundf(self)
        }
    }

    #[inline]
    fn trunc(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.trunc()
        }
        #[cfg(not(feature = "std"))]
        {
            truncf(self)
        }
    }

    #[inline]
    fn fract(self) -> f32 {
        #[cfg(feature = "std")]
        {
            self.fract()
        }
        // libm provides no fract: inherent is defined as self - self.trunc(),
        // reproduced exactly (including NaN, infinite, and signed-zero edges).
        #[cfg(not(feature = "std"))]
        {
            self - truncf(self)
        }
    }

    #[inline]
    fn gamma(self) -> f32 {
        SpecialFunc::gamma(self)
    }

    #[inline]
    fn lgamma(self) -> f32 {
        SpecialFunc::lgamma(self)
    }

    #[inline]
    fn digamma(self) -> f32 {
        SpecialFunc::digamma(self)
    }

    #[inline]
    fn trigamma(self) -> f32 {
        SpecialFunc::trigamma(self)
    }

    #[inline]
    fn tetragamma(self) -> f32 {
        SpecialFunc::tetragamma(self)
    }

    #[inline]
    fn polygamma(self, n: i32) -> f32 {
        SpecialFunc::polygamma_n(n, self)
    }

    #[inline]
    fn erf(self) -> f32 {
        SpecialFunc::erf(self)
    }

    #[inline]
    fn erfc(self) -> f32 {
        SpecialFunc::erfc(self)
    }

    #[inline]
    fn zeta(self) -> f32 {
        SpecialFunc::zeta(self)
    }

    #[inline]
    fn zeta_deriv(self, n: i32) -> f32 {
        SpecialFunc::zeta_deriv(n, self)
    }

    #[inline]
    fn bessel_j(self, n: i32) -> f32 {
        SpecialFunc::besselj(n, self)
    }

    #[inline]
    fn bessel_y(self, n: i32) -> f32 {
        SpecialFunc::bessely(n, self)
    }

    #[inline]
    fn bessel_i(self, n: i32) -> f32 {
        SpecialFunc::besseli(n, self)
    }

    #[inline]
    fn bessel_k(self, n: i32) -> f32 {
        SpecialFunc::besselk(n, self)
    }

    #[inline]
    fn beta(self, b: f32) -> f32 {
        SpecialFunc::beta(self, b)
    }

    #[inline]
    fn elliptic_k(self) -> f32 {
        SpecialFunc::elliptic_k(self)
    }

    #[inline]
    fn elliptic_e(self) -> f32 {
        SpecialFunc::elliptic_e(self)
    }

    #[inline]
    fn hermite(self, n: i32) -> f32 {
        SpecialFunc::hermite(n, self)
    }

    #[inline]
    fn lambert_w0(self) -> f32 {
        SpecialFunc::lambertw0(self)
    }

    #[inline]
    fn lambert_wm1(self) -> f32 {
        SpecialFunc::lambertwm1(self)
    }

    #[inline]
    fn assoc_legendre(self, l: i32, m: i32) -> f32 {
        SpecialFunc::assoc_legendre(l, m, self)
    }

    #[inline]
    fn spherical_harmonic(self, l: i32, m: i32, phi: f32) -> f32 {
        SpecialFunc::spherical_harmonic(l, m, self, phi)
    }
}

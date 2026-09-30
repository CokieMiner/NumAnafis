//! Special mathematical functions — generic over float and integer types.
//!
//! Each sub-module works via the `SpecFloat` and `SpecInt` traits so that
//! algorithms are backend-switchable with precision-specific constants
//! behind compilation walls.
//!
//! # Bibliography
//!
//! References cited throughout the `special` sub-modules follow this key:
//! - **\[Abramowitz64\]** Abramowitz, M. & Stegun, I.A. (1964). *Handbook of Mathematical Functions*. National Bureau of Standards.
//! - **\[Borwein87\]** Borwein, J. M., & Borwein, P. B. (1987). *Pi and the AGM*. Wiley.
//! - **\[Borwein00\]** Borwein, J. M., Bradley, D. M., & Crandall, R. E. (2000). "Computational strategies for the Riemann zeta function." *J. Comput. Appl. Math.*, 121(1-2), 247-296.
//! - **\[Carlson95\]** Carlson, B. C. (1995). "Numerical computation of real or complex elliptic integrals." *Numerical Algorithms*, 10(1), 13-26.
//! - **\[Cephes\]** Moshier, S. L. (1989). *Methods and Programs for Mathematical Functions*. Ellis Horwood Limited. (Cephes Mathematical Library)
//! - **\[Clenshaw55\]** Clenshaw, C. W. (1955). "A note on the summation of Chebyshev series." *Math. Tables Aids Comput.*, 9(51), 118-120.
//! - **\[Corless96\]** Corless, R. M., et al. (1996). "On the Lambert W function." *Adv. Comput. Math.*, 5(1), 329-359.
//! - **\[DLMF\]** NIST Digital Library of Mathematical Functions. `<https://dlmf.nist.gov/>`, Release 1.1.12 of 2023-12-15. F. W. J. Olver et al., eds.
//! - **\[GKP94\]** Graham, R. L., Knuth, D. E., & Patashnik, O. (1994). *Concrete Mathematics* (2nd ed.). Addison-Wesley.
//! - **\[Lanczos64\]** Lanczos, C. (1964). "A precision approximation of the gamma function." *J. SIAM Numer. Anal. Ser. B*, 1(1), 86-96.
//! - **\[Miller52\]** Miller, J. C. P. (1952). "A method for the determination of converging factors..." *Proc. Camb. Philos. Soc.*, 48(2), 243-254.
mod bessel_ik;
mod bessel_jy;
mod beta;
mod coeffs_f32;
mod coeffs_f64;
mod complex;
mod elliptic;
mod erf;
mod funcs;
mod gamma;
mod helpers;
mod hermite;
mod impls;
mod lambert_w;
mod legendre;
mod polygamma;
mod traits;
mod zeta;
mod zeta_deriv;

pub use complex::{C, ComplexMath};
pub use funcs::SpecialFunc;
pub use helpers::SpecialHelper;
pub use traits::{SpecFloat, SpecInt};

//! The state the model's pieces share while one evaluation runs.
//!
//! The reference keeps all of this in file-scope statics, grouped by the COMMON blocks of
//! the Fortran it was translated from: PARMB, DMIX, MESO7 and LPOLY. Here it is one struct,
//! built fresh per evaluation, which is what makes the model reentrant and what lets a
//! caller run one per thread.
//!
//! # Where that changes an answer
//!
//! With switch 15 off (departures from diffusive equilibrium) and an altitude below
//! 72.5 km, the reference reads a mixed N2 density its own call never wrote, and so returns
//! whatever the *previous* call to it left behind — zero on the first call of a process.
//! Starting from zero is the behavior this reproduces. Nothing else in the model depends on
//! what a previous call did.

use crate::scalar::Scalar;

/// Everything one evaluation carries between its pieces.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Work<T: Scalar> {
    /// Gravity at the surface at this latitude, cm/s^2. `gsurf` in the reference.
    pub gsurf: T,
    /// The effective Earth radius at this latitude, km. `re` in the reference.
    pub re: T,

    /// Mixed densities at the altitude asked for: He, O, N2, O2, Ar, H and N.
    ///
    /// `dm04` and friends, the DMIX block. Only N2 is read outside the species that set
    /// it, and only below 72.5 km.
    pub dm04: T,
    pub dm16: T,
    pub dm28: T,
    pub dm32: T,
    pub dm40: T,
    pub dm01: T,
    pub dm14: T,

    /// Lower thermosphere temperature nodes, 72.5 to 120 km. `meso_tn1`.
    pub tn1: [T; 5],
    /// Upper stratosphere and mesosphere nodes, 32.5 to 72.5 km. `meso_tn2`.
    pub tn2: [T; 4],
    /// Troposphere and lower stratosphere nodes, 0 to 32.5 km. `meso_tn3`.
    pub tn3: [T; 5],
    /// Temperature gradients at the ends of the `tn1` spline. `meso_tgn1`.
    pub tgn1: [T; 2],
    /// Temperature gradients at the ends of the `tn2` spline. `meso_tgn2`.
    pub tgn2: [T; 2],
    /// Temperature gradients at the ends of the `tn3` spline. `meso_tgn3`.
    pub tgn3: [T; 2],

    /// Whether [`crate::expansion::prepare`] has run for the input in hand.
    ///
    /// Only read by a debug assertion. The terms it computes — the Legendre polynomials,
    /// the local time harmonics and the flux departure — are functions of the input alone,
    /// so they are the same on every one of `globe7`'s eleven calls per evaluation, and
    /// computing them once is the single largest saving in the model. The cost is that a
    /// `Work` now belongs to one input, where the reference recomputed them every call and
    /// so belonged to none.
    pub prepared: bool,
    /// F10.7 average less 150. `dfa`.
    pub dfa: T,
    /// Associated Legendre polynomials of the latitude, `plg[m][n]`.
    pub plg: [[T; 9]; 4],
    /// Sine and cosine of the day of year at one and two cycles a year.
    ///
    /// Only the `fast` path reads these, and they are zero-sized without it. The
    /// expansions ask for the cosine of the day of year at eight different phase shifts
    /// between them, and from these two each is an angle-sum away.
    pub doy1: crate::expansion::Harmonic<f64>,
    pub doy2: crate::expansion::Harmonic<f64>,
    /// Sine and cosine of the local solar time, universal time and longitude at the rates
    /// the phase-shifted terms use. Zero-sized without `fast`, as above.
    pub hour: crate::expansion::Harmonic<T>,
    pub universal: crate::expansion::Harmonic<T>,
    pub east: crate::expansion::Harmonic<T>,
    pub east2: crate::expansion::Harmonic<T>,
    /// Sine and cosine of the local solar time, at one, two and three cycles a day.
    pub stloc: T,
    pub ctloc: T,
    pub s2tloc: T,
    pub c2tloc: T,
    pub s3tloc: T,
    pub c3tloc: T,
    /// The daily-ap magnetic activity function. `apdf`.
    pub apdf: T,
    /// The last set of three-hour magnetic activity functions computed, and the two
    /// coefficients they belong to.
    ///
    /// `g0` is a function of two coefficients and one ap value, and nine of the eleven
    /// tables `globe7` is called with carry the same two. Keeping the last set means the
    /// six exponentials it needs are computed three times per evaluation rather than ten,
    /// which is the largest saving available without changing an answer. The values are
    /// the same ones the reference computes, from the same expression.
    pub magnetic_key: [f64; 2],
    pub magnetic_cached: bool,
    pub magnetic_g0: [T; 6],
    /// The ap-history magnetic activity functions. `apt`; only the first is ever set.
    pub apt: [T; 4],
}

impl<T: Scalar> Work<T> {
    /// All zero, which is where the reference's statics start.
    pub fn new() -> Self {
        let zero = T::lift(0.0);
        Self {
            gsurf: zero,
            re: zero,
            dm04: zero,
            dm16: zero,
            dm28: zero,
            dm32: zero,
            dm40: zero,
            dm01: zero,
            dm14: zero,
            tn1: [zero; 5],
            tn2: [zero; 4],
            tn3: [zero; 5],
            tgn1: [zero; 2],
            tgn2: [zero; 2],
            tgn3: [zero; 2],
            prepared: false,
            dfa: zero,
            plg: [[zero; 9]; 4],
            doy1: crate::expansion::Harmonic::none(),
            doy2: crate::expansion::Harmonic::none(),
            hour: crate::expansion::Harmonic::none(),
            universal: crate::expansion::Harmonic::none(),
            east: crate::expansion::Harmonic::none(),
            east2: crate::expansion::Harmonic::none(),
            stloc: zero,
            ctloc: zero,
            s2tloc: zero,
            c2tloc: zero,
            s3tloc: zero,
            c3tloc: zero,
            apdf: zero,
            magnetic_key: [f64::NAN; 2],
            magnetic_cached: false,
            magnetic_g0: [zero; 6],
            apt: [zero; 4],
        }
    }
}

//! The expansions in latitude, day of year, local solar time, longitude, flux and magnetic
//! activity: `globe7` for the thermosphere and `glob7s` for below it.
//!
//! These are where the model's three thousand coefficients are spent. Each builds a sum of
//! fourteen terms, one per family of variation, and each family is scaled by the switch
//! that governs it. The switches are indexed here by the reference's own numbers; the
//! names are on [`crate::Switches`].
//!
//! `globe7` leaves the Legendre polynomials of the latitude, the local time harmonics, the
//! flux departure and the magnetic activity functions in [`Work`], and `glob7s` reads them
//! rather than recomputing them. So `globe7` must be called first, which the model does.
//!
//! Terms that are products of a day-of-year cosine and a coefficient stay `f64`, because
//! the day of year is an integer: `cd32`, `cd18`, `cd14` and `cd39` are real constants of
//! the call and cost nothing to carry at a dual type.

use crate::helpers::{DGTR, DR, HR, SR};
use crate::input::{Input, MagneticActivity};
use crate::scalar::Scalar;
use crate::switches::Resolved;
use crate::work::Work;

/// How many phase-shifted cosines a coefficient table supplies.
///
/// Every cosine in the expansions but one has the form `cos(w * (x - p))`: a rate `w` and
/// a quantity `x` that the whole evaluation shares, against a phase `p` that is one of the
/// table's own coefficients. There are eighteen such coefficients between the two
/// expansions, across six rates.
pub(crate) const PHASES: usize = 18;

/// The coefficient index and angular rate of each, in the order [`Phase`] arrays hold them.
///
/// The first eight are day-of-year phases, which both expansions read; the rest are read
/// only by [`globe7`], and four of them sit past the end of the shorter tables [`glob7s`]
/// is called with, so [`Phase::of`] leaves those empty.
const PHASE_OF: [(usize, f64); PHASES] = [
    (13, DR),
    (31, DR),
    (81, DR),
    (84, DR),
    (17, 2.0 * DR),
    (38, 2.0 * DR),
    (86, 2.0 * DR),
    (88, 2.0 * DR),
    (131, HR),
    (124, HR),
    (71, SR),
    (79, SR),
    (58, SR),
    (75, SR),
    (97, DGTR),
    (136, DGTR),
    (63, DGTR),
    (118, DGTR),
];

/// Where each coefficient sits in a [`Phase`] array.
pub(crate) const P13: usize = 0;
pub(crate) const P31: usize = 1;
pub(crate) const P81: usize = 2;
pub(crate) const P84: usize = 3;
pub(crate) const P17: usize = 4;
pub(crate) const P38: usize = 5;
pub(crate) const P86: usize = 6;
pub(crate) const P88: usize = 7;
const P131: usize = 8;
const P124: usize = 9;
const P71: usize = 10;
const P79: usize = 11;
const P58: usize = 12;
const P75: usize = 13;
const P97: usize = 14;
const P136: usize = 15;
const P63: usize = 16;
const P118: usize = 17;

/// A phase of zero, for a cosine the reference writes without a shift.
const UNSHIFTED: Phase = Phase::zero();

/// The sine and cosine of one phase shift, which the `fast` path needs and the exact path
/// does not.
///
/// Zero-sized when the feature is off, so neither the model nor the arithmetic carries
/// anything for a path that is not taken.
#[cfg(feature = "fast")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Phase {
    sin: f64,
    cos: f64,
}

/// The sine and cosine of one phase shift. See the other definition.
#[cfg(not(feature = "fast"))]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Phase;

impl Phase {
    /// The phases of one coefficient table, worked out once when a model is built.
    ///
    /// A table shorter than an index asked for does not use that phase, and gets the
    /// unshifted one: the shorter tables are the hundred-coefficient rows `glob7s` reads,
    /// and the phases past that point belong to terms only `globe7` has.
    pub(crate) fn of(p: &[f64]) -> [Self; PHASES] {
        let mut phases = [Self::zero(); PHASES];
        for (phase, (index, omega)) in phases.iter_mut().zip(PHASE_OF) {
            if index < p.len() {
                *phase = Self::at(omega * p[index]);
            }
        }
        phases
    }

    #[cfg(feature = "fast")]
    fn at(angle: f64) -> Self {
        Self {
            sin: angle.sin(),
            cos: angle.cos(),
        }
    }

    #[cfg(not(feature = "fast"))]
    fn at(_angle: f64) -> Self {
        Self
    }

    /// The phase of an unshifted cosine, whose sine is zero and cosine one.
    #[cfg(feature = "fast")]
    const fn zero() -> Self {
        Self { sin: 0.0, cos: 1.0 }
    }

    /// The phase of an unshifted cosine. See the other definition.
    #[cfg(not(feature = "fast"))]
    const fn zero() -> Self {
        Self
    }
}

/// The sine and cosine of an angle the whole evaluation shares.
///
/// Zero-sized when `fast` is off, like [`Phase`].
#[cfg(feature = "fast")]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Harmonic<T> {
    sin: T,
    cos: T,
}

/// The sine and cosine of an angle. See the other definition.
#[cfg(not(feature = "fast"))]
#[derive(Clone, Copy, Debug)]
pub(crate) struct Harmonic<T>(core::marker::PhantomData<T>);

impl<T: Scalar> Harmonic<T> {
    #[cfg(feature = "fast")]
    pub(crate) fn of(angle: T) -> Self {
        Self {
            sin: angle.sin(),
            cos: angle.cos(),
        }
    }

    #[cfg(not(feature = "fast"))]
    pub(crate) fn of(_angle: T) -> Self {
        Self(core::marker::PhantomData)
    }

    #[cfg(not(feature = "fast"))]
    pub(crate) fn none() -> Self {
        Self(core::marker::PhantomData)
    }

    #[cfg(feature = "fast")]
    pub(crate) fn none() -> Self {
        Self {
            sin: T::lift(0.0),
            cos: T::lift(0.0),
        }
    }
}

/// `cos(omega * (x - p))`.
///
/// Exactly, by calling `cos`, or under `fast` by the angle-sum identity from the
/// evaluation's own `cos(omega * x)` and `sin(omega * x)` and the phase's own sine and
/// cosine. The identity is exact in real arithmetic and is two roundings rather than one
/// in floating point, so the two differ in the last bits; what it buys is one `cos` call
/// per evaluation where the expansions ask for a hundred.
#[inline(always)]
#[allow(unused_variables)]
fn shifted<T: Scalar>(angle: Harmonic<T>, phase: Phase, omega: f64, x: T, p: f64) -> T {
    #[cfg(feature = "fast")]
    return angle.cos * phase.cos + angle.sin * phase.sin;
    #[cfg(not(feature = "fast"))]
    return ((x - p) * omega).cos();
}

/// `cos(omega * (x - p) + two_omega * y)`, the one cosine of a sum the expansions ask for.
///
/// Expanded twice: once to split the sum, and once for each phase-shifted half.
#[inline(always)]
#[allow(unused_variables, clippy::too_many_arguments)]
fn mixed<T: Scalar>(
    angle: Harmonic<T>,
    doubled: Harmonic<T>,
    phase: Phase,
    x: T,
    p: f64,
    y: T,
) -> T {
    #[cfg(feature = "fast")]
    return shifted(angle, phase, SR, x, p) * doubled.cos
        - shifted_sin(angle, phase, SR, x, p) * doubled.sin;
    #[cfg(not(feature = "fast"))]
    return ((x - p) * SR + y * (2.0 * DGTR)).cos();
}

/// `sin(omega * (x - p))`, the counterpart of [`shifted`].
#[inline(always)]
#[allow(unused_variables)]
fn shifted_sin<T: Scalar>(angle: Harmonic<T>, phase: Phase, omega: f64, x: T, p: f64) -> T {
    #[cfg(feature = "fast")]
    return angle.sin * phase.cos - angle.cos * phase.sin;
    #[cfg(not(feature = "fast"))]
    return ((x - p) * omega).sin();
}

/// The three-hour magnetic activity function, equation A24d.
fn g0<T: Scalar>(a: T, p: &[f64]) -> T {
    // sqrt(p[24]*p[24]) is the reference's way of writing |p[24]|.
    let p24 = p[24].abs();
    (a - 4.0) + ((a - 4.0) + (((a - 4.0) * (-p24)).exp() - 1.0) / p24) * (p[25] - 1.0)
}

/// `g0` of the six three-hour ap values, from the cache when the coefficients match.
///
/// It is a function of `p[24]`, `p[25]` and the ap history, and the history does not
/// change within an evaluation, so the only thing that can invalidate the cache is a
/// coefficient table with different values in those two slots. Two of the eleven have
/// them, so this computes six exponentials three times rather than ten.
fn g0_history<'a, T: Scalar>(p: &[f64], ap: &[T; 7], work: &'a mut Work<T>) -> &'a [T; 6] {
    let key = [p[24], p[25]];
    if !work.magnetic_cached || work.magnetic_key != key {
        for (slot, value) in work.magnetic_g0.iter_mut().enumerate() {
            *value = g0(ap[slot + 1], p);
        }
        work.magnetic_key = key;
        work.magnetic_cached = true;
    }
    &work.magnetic_g0
}

/// The normalizing sum of equation A24c.
fn sumex<T: Scalar>(ex: T) -> T {
    (T::lift(1.0) - powers(ex).nineteen) / (T::lift(1.0) - ex) * ex.sqrt() + 1.0
}

/// The powers of `ex` that the magnetic activity functions need.
///
/// The reference asks for the 3rd, 4th, 8th, 12th and 19th powers, each as a `pow` call
/// with a constant exponent. Built by repeated squaring instead they cost seven
/// multiplications between them rather than five calls into the maths library, which is
/// worth about a sixth of an evaluation in the ap-history mode and nothing in the other. It is not the same arithmetic: `pow` is
/// correctly rounded and a product of products is not, so the last bits differ. That is
/// what the `fast` feature chooses.
struct Powers<T> {
    three: T,
    four: T,
    eight: T,
    twelve: T,
    nineteen: T,
}

#[cfg(feature = "fast")]
fn powers<T: Scalar>(ex: T) -> Powers<T> {
    let two = ex * ex;
    let three = two * ex;
    let four = two * two;
    let eight = four * four;
    let sixteen = eight * eight;
    Powers {
        three,
        four,
        eight,
        twelve: eight * four,
        nineteen: sixteen * three,
    }
}

#[cfg(not(feature = "fast"))]
fn powers<T: Scalar>(ex: T) -> Powers<T> {
    Powers {
        three: ex.powf(3.0),
        four: ex.powf(4.0),
        eight: ex.powf(8.0),
        twelve: ex.powf(12.0),
        nineteen: ex.powf(19.0),
    }
}

/// The weighted history of magnetic activity, equation A24a.
///
/// `g` holds `g0` of `ap[1]` through `ap[6]`, in that order.
fn sg0<T: Scalar>(ex: T, g: &[T; 6]) -> T {
    let power = powers(ex);
    (g[0]
        + (g[1] * ex
            + g[2] * ex * ex
            + g[3] * power.three
            + (g[4] * power.four + g[5] * power.twelve) * (T::lift(1.0) - power.eight)
                / (T::lift(1.0) - ex)))
        / sumex(ex)
}

/// Computes the parts of the expansions that depend only on the input.
///
/// The Legendre polynomials of the latitude, the six local time harmonics and the flux
/// departure are the same for every coefficient table, so they are computed once per
/// evaluation here rather than on each of [`globe7`]'s eleven calls. That is worth about a
/// seventh of an evaluation, and the values are the reference's own.
///
/// A [`Work`] therefore belongs to one input. [`globe7`] and [`glob7s`] assert in debug
/// builds that this has run.
pub(crate) fn prepare<T: Scalar>(input: &Input<T>, switches: &Resolved, work: &mut Work<T>) {
    let sw = &switches.sw;
    let tloc = input.local_solar_time_h;

    // The associated Legendre polynomials of the latitude. The reference leaves three of
    // them commented out, and those stay out: plg[0][7], plg[1][7] and plg[1][8].
    let c = (input.latitude_deg * DGTR).sin();
    let s = (input.latitude_deg * DGTR).cos();
    let c2 = c * c;
    let c4 = c2 * c2;
    let s2 = s * s;

    let plg = &mut work.plg;
    plg[0][1] = c;
    plg[0][2] = (c2 * 3.0 - 1.0) * 0.5;
    plg[0][3] = (c * 5.0 * c2 - c * 3.0) * 0.5;
    plg[0][4] = (c4 * 35.0 - c2 * 30.0 + 3.0) / 8.0;
    plg[0][5] = (c2 * 63.0 * c2 * c - c2 * 70.0 * c + c * 15.0) / 8.0;
    plg[0][6] = (c * 11.0 * plg[0][5] - plg[0][4] * 5.0) / 6.0;
    plg[1][1] = s;
    plg[1][2] = c * 3.0 * s;
    plg[1][3] = (c2 * 5.0 - 1.0) * 1.5 * s;
    plg[1][4] = (c2 * 7.0 * c - c * 3.0) * 2.5 * s;
    plg[1][5] = (c4 * 21.0 - c2 * 14.0 + 1.0) * 1.875 * s;
    plg[1][6] = (c * 11.0 * plg[1][5] - plg[1][4] * 6.0) / 5.0;
    plg[2][2] = s2 * 3.0;
    plg[2][3] = s2 * 15.0 * c;
    plg[2][4] = (c2 * 7.0 - 1.0) * 7.5 * s2;
    plg[2][5] = c * 3.0 * plg[2][4] - plg[2][3] * 2.0;
    plg[2][6] = (c * 11.0 * plg[2][5] - plg[2][4] * 7.0) / 4.0;
    plg[2][7] = (c * 13.0 * plg[2][6] - plg[2][5] * 8.0) / 5.0;
    plg[3][3] = s2 * 15.0 * s;
    plg[3][4] = s2 * 105.0 * s * c;
    plg[3][5] = (c * 9.0 * plg[3][4] - plg[3][3] * 7.0) / 2.0;
    plg[3][6] = (c * 11.0 * plg[3][5] - plg[3][4] * 8.0) / 3.0;

    if !(sw[7] == 0.0 && sw[8] == 0.0 && sw[14] == 0.0) {
        work.stloc = (tloc * HR).sin();
        work.ctloc = (tloc * HR).cos();
        work.s2tloc = (tloc * (2.0 * HR)).sin();
        work.c2tloc = (tloc * (2.0 * HR)).cos();
        work.s3tloc = (tloc * (3.0 * HR)).sin();
        work.c3tloc = (tloc * (3.0 * HR)).cos();
    }

    // The day of year at one and two cycles a year, which the `fast` path turns every
    // phase-shifted day-of-year cosine into an angle-sum from. Nothing without it.
    let doy = f64::from(input.day_of_year);
    work.doy1 = Harmonic::of(DR * doy);
    work.doy2 = Harmonic::of(2.0 * DR * doy);
    work.hour = Harmonic::of(input.local_solar_time_h * HR);
    work.universal = Harmonic::of(input.seconds_of_day * SR);
    work.east = Harmonic::of(input.longitude_deg * DGTR);
    work.east2 = Harmonic::of(input.longitude_deg * (2.0 * DGTR));

    // The flux departure belongs here too: glob7s reads it, and it depends only on the
    // input.
    work.dfa = input.f107_average - 150.0;

    // The magnetic activity functions depend on the ap values, which belong to this input,
    // so anything cached from a previous one is stale.
    work.magnetic_cached = false;
    work.prepared = true;
}

/// The expansion for the upper thermosphere: `globe7`.
///
/// Returns the sum of the fourteen term families for the coefficient set `p`, which is one
/// of `PT`, `PS` or a row of `PD`.
pub(crate) fn globe7<T: Scalar>(
    p: &[f64],
    phase: &[Phase; PHASES],
    input: &Input<T>,
    switches: &Resolved,
    work: &mut Work<T>,
) -> T {
    debug_assert!(work.prepared, "prepare() has not run for this input");
    let (sw, swc) = (&switches.sw, &switches.swc);
    let mut t = [T::lift(0.0); 14];
    let tloc = input.local_solar_time_h;
    let doy = f64::from(input.day_of_year);

    // Day-of-year cosines. Real constants of the call, since the day is an integer.
    let cd32 = shifted(work.doy1, phase[P31], DR, doy, p[31]);
    let cd18 = shifted(work.doy2, phase[P17], 2.0 * DR, doy, p[17]);
    let cd14 = shifted(work.doy1, phase[P13], DR, doy, p[13]);
    let cd39 = shifted(work.doy2, phase[P38], 2.0 * DR, doy, p[38]);

    let df = input.f107 - input.f107_average;
    let dfa = work.dfa;
    t[0] = df * p[19] * (dfa * p[59] + 1.0) + df * p[20] * df + dfa * p[21] + dfa * dfa * p[29];
    let f1 = (dfa * p[47] + df * p[19] + df * p[20] * df) * swc[1] + 1.0;
    let f2 = (dfa * p[49] + df * p[19] + df * p[20] * df) * swc[1] + 1.0;

    // Copied rather than read through the struct on every term: the local array lets the
    // values stay in registers, and measured it is a percent faster that way.
    let plg = work.plg;

    // Time independent.
    t[1] = plg[0][2] * p[1]
        + plg[0][4] * p[2]
        + plg[0][6] * p[22]
        + plg[0][2] * p[14] * dfa * swc[1]
        + plg[0][1] * p[26];

    // Symmetrical annual.
    t[2] = T::lift(p[18] * cd32);

    // Symmetrical semiannual.
    t[3] = (plg[0][2] * p[16] + p[15]) * cd18;

    // Asymmetrical annual.
    t[4] = f1 * (plg[0][1] * p[9] + plg[0][3] * p[10]) * cd14;

    // Asymmetrical semiannual.
    t[5] = plg[0][1] * p[37] * cd39;

    // Diurnal.
    if sw[7] != 0.0 {
        let t71 = plg[1][2] * p[11] * cd14 * swc[5];
        let t72 = plg[1][2] * p[12] * cd14 * swc[5];
        t[6] = f2
            * ((plg[1][1] * p[3] + plg[1][3] * p[4] + plg[1][5] * p[27] + t71) * work.ctloc
                + (plg[1][1] * p[6] + plg[1][3] * p[7] + plg[1][5] * p[28] + t72) * work.stloc);
    }

    // Semidiurnal.
    if sw[8] != 0.0 {
        let t81 = (plg[2][3] * p[23] + plg[2][5] * p[35]) * cd14 * swc[5];
        let t82 = (plg[2][3] * p[33] + plg[2][5] * p[36]) * cd14 * swc[5];
        t[7] = f2
            * ((plg[2][2] * p[5] + plg[2][4] * p[41] + t81) * work.c2tloc
                + (plg[2][2] * p[8] + plg[2][4] * p[42] + t82) * work.s2tloc);
    }

    // Terdiurnal.
    if sw[14] != 0.0 {
        t[13] = f2
            * ((plg[3][3] * p[39] + (plg[3][4] * p[93] + plg[3][6] * p[46]) * cd14 * swc[5])
                * work.s3tloc
                + (plg[3][3] * p[40] + (plg[3][4] * p[94] + plg[3][6] * p[48]) * cd14 * swc[5])
                    * work.c3tloc);
    }

    // Magnetic activity, in whichever of its two modes the input selects.
    let sw9 = input.magnetic.sw();
    match &input.magnetic {
        MagneticActivity::History { ap } => {
            if p[51] != 0.0 {
                let mut exp1 = (T::lift(-10800.0 * p[51].abs())
                    / ((T::lift(45.0) - input.latitude_deg.abs()) * p[138] + 1.0))
                    .exp();
                if exp1.re() > 0.99999 {
                    exp1 = T::lift(0.99999);
                }
                // The reference clamps p[24] up to 1e-4 here, writing into the coefficient
                // table. No table it is ever called with has a value that small, so the
                // clamp cannot fire; `test_the_p24_clamp_cannot_fire` holds that down.
                let g = *g0_history(p, ap, work);
                work.apt[0] = sg0(exp1, &g);
                if sw9 != 0.0 {
                    t[8] = work.apt[0]
                        * (T::lift(p[50])
                            + plg[0][2] * p[96]
                            + plg[0][4] * p[54]
                            + (plg[0][1] * p[125] + plg[0][3] * p[126] + plg[0][5] * p[127])
                                * cd14
                                * swc[5]
                            + (plg[1][1] * p[128] + plg[1][3] * p[129] + plg[1][5] * p[130])
                                * swc[7]
                                * shifted(work.hour, phase[P131], HR, tloc, p[131]));
                }
            }
        }
        _ => {
            let apd = input.magnetic.daily() - 4.0;
            let p44 = if p[43] < 0.0 { 1.0e-5 } else { p[43] };
            let p45 = p[44];
            work.apdf = apd + (apd + ((apd * (-p44)).exp() - 1.0) / p44) * (p45 - 1.0);
            if sw9 != 0.0 {
                t[8] = work.apdf
                    * (T::lift(p[32])
                        + plg[0][2] * p[45]
                        + plg[0][4] * p[34]
                        + (plg[0][1] * p[100] + plg[0][3] * p[101] + plg[0][5] * p[102])
                            * cd14
                            * swc[5]
                        + (plg[1][1] * p[121] + plg[1][3] * p[122] + plg[1][5] * p[123])
                            * swc[7]
                            * shifted(work.hour, phase[P124], HR, tloc, p[124]));
            }
        }
    }

    if sw[10] != 0.0 && input.longitude_deg.re() > -1000.0 {
        // Longitudinal.
        if sw[11] != 0.0 {
            t[10] = (dfa * p[80] * swc[1] + 1.0)
                * ((plg[1][2] * p[64]
                    + plg[1][4] * p[65]
                    + plg[1][6] * p[66]
                    + plg[1][1] * p[103]
                    + plg[1][3] * p[104]
                    + plg[1][5] * p[105]
                    + (plg[1][1] * p[109] + plg[1][3] * p[110] + plg[1][5] * p[111])
                        * swc[5]
                        * cd14)
                    * shifted(work.east, UNSHIFTED, DGTR, input.longitude_deg, 0.0)
                    + (plg[1][2] * p[90]
                        + plg[1][4] * p[91]
                        + plg[1][6] * p[92]
                        + plg[1][1] * p[106]
                        + plg[1][3] * p[107]
                        + plg[1][5] * p[108]
                        + (plg[1][1] * p[112] + plg[1][3] * p[113] + plg[1][5] * p[114])
                            * swc[5]
                            * cd14)
                        * shifted_sin(work.east, UNSHIFTED, DGTR, input.longitude_deg, 0.0));
        }

        // Universal time, and mixed universal time and longitude.
        if sw[12] != 0.0 {
            t[11] = (plg[0][1] * p[95] + 1.0)
                * (dfa * p[81] * swc[1] + 1.0)
                * (plg[0][1] * p[119] * swc[5] * cd14 + 1.0)
                * ((plg[0][1] * p[68] + plg[0][3] * p[69] + plg[0][5] * p[70])
                    * shifted(work.universal, phase[P71], SR, input.seconds_of_day, p[71]));
            t[11] += (plg[2][3] * p[76] + plg[2][5] * p[77] + plg[2][7] * p[78])
                * swc[11]
                * mixed(
                    work.universal,
                    work.east2,
                    phase[P79],
                    input.seconds_of_day,
                    p[79],
                    input.longitude_deg,
                )
                * (dfa * p[137] * swc[1] + 1.0);
        }

        // Mixed magnetic activity, universal time and longitude.
        if sw[13] != 0.0 {
            match &input.magnetic {
                MagneticActivity::History { .. } => {
                    if p[51] != 0.0 {
                        t[12] = work.apt[0]
                            * swc[11]
                            * (plg[0][1] * p[132] + 1.0)
                            * ((plg[1][2] * p[52] + plg[1][4] * p[98] + plg[1][6] * p[67])
                                * shifted(work.east, phase[P97], DGTR, input.longitude_deg, p[97]))
                            + work.apt[0]
                                * swc[11]
                                * swc[5]
                                * (plg[1][1] * p[133] + plg[1][3] * p[134] + plg[1][5] * p[135])
                                * cd14
                                * shifted(
                                    work.east,
                                    phase[P136],
                                    DGTR,
                                    input.longitude_deg,
                                    p[136],
                                )
                            + work.apt[0]
                                * swc[12]
                                * (plg[0][1] * p[55] + plg[0][3] * p[56] + plg[0][5] * p[57])
                                * shifted(
                                    work.universal,
                                    phase[P58],
                                    SR,
                                    input.seconds_of_day,
                                    p[58],
                                );
                    }
                }
                _ => {
                    t[12] = work.apdf
                        * swc[11]
                        * (plg[0][1] * p[120] + 1.0)
                        * ((plg[1][2] * p[60] + plg[1][4] * p[61] + plg[1][6] * p[62])
                            * shifted(work.east, phase[P63], DGTR, input.longitude_deg, p[63]))
                        + work.apdf
                            * swc[11]
                            * swc[5]
                            * (plg[1][1] * p[115] + plg[1][3] * p[116] + plg[1][5] * p[117])
                            * cd14
                            * shifted(work.east, phase[P118], DGTR, input.longitude_deg, p[118])
                        + work.apdf
                            * swc[12]
                            * (plg[0][1] * p[83] + plg[0][3] * p[84] + plg[0][5] * p[85])
                            * shifted(work.universal, phase[P75], SR, input.seconds_of_day, p[75]);
                }
            }
        }
    }

    // Coefficients 82, 89, 99 and 139 to 149 are not used.
    let mut tinf = T::lift(p[30]);
    for (index, term) in t.iter().enumerate() {
        let switch = if index == 8 { sw9 } else { sw[index + 1] };
        tinf += *term * switch.abs();
    }
    tinf
}

/// The expansion for the lower atmosphere: `glob7s`.
///
/// Reads the state `globe7` left in `work` rather than recomputing it, which is why the
/// model calls that first. `p` is a row of `PTL` or of `PMA`.
pub(crate) fn glob7s<T: Scalar>(
    p: &[f64],
    phase: &[Phase; PHASES],
    input: &Input<T>,
    switches: &Resolved,
    work: &Work<T>,
) -> T {
    let (sw, swc) = (&switches.sw, &switches.swc);
    let mut t = [T::lift(0.0); 14];
    let doy = f64::from(input.day_of_year);
    let plg = &work.plg;

    // The reference checks that p[99] is its parameter-set marker of 2, writing it in if
    // the slot is zero. Every row of PTL and PMA already carries 2, so neither the write
    // nor the complaint can happen; `test_every_lower_atmosphere_table_is_parameter_set_2`
    // holds that down.
    debug_assert_eq!(p[99], 2.0, "glob7s called with the wrong parameter set");

    let cd32 = shifted(work.doy1, phase[P31], DR, doy, p[31]);
    let cd18 = shifted(work.doy2, phase[P17], 2.0 * DR, doy, p[17]);
    let cd14 = shifted(work.doy1, phase[P13], DR, doy, p[13]);
    let cd39 = shifted(work.doy2, phase[P38], 2.0 * DR, doy, p[38]);

    // F10.7.
    t[0] = work.dfa * p[21];

    // Time independent.
    t[1] = plg[0][2] * p[1]
        + plg[0][4] * p[2]
        + plg[0][6] * p[22]
        + plg[0][1] * p[26]
        + plg[0][3] * p[14]
        + plg[0][5] * p[59];

    // Symmetrical annual.
    t[2] = (T::lift(p[18]) + plg[0][2] * p[47] + plg[0][4] * p[29]) * cd32;

    // Symmetrical semiannual.
    t[3] = (T::lift(p[15]) + plg[0][2] * p[16] + plg[0][4] * p[30]) * cd18;

    // Asymmetrical annual.
    t[4] = (plg[0][1] * p[9] + plg[0][3] * p[10] + plg[0][5] * p[20]) * cd14;

    // Asymmetrical semiannual.
    t[5] = plg[0][1] * p[37] * cd39;

    // Diurnal.
    if sw[7] != 0.0 {
        let t71 = plg[1][2] * p[11] * cd14 * swc[5];
        let t72 = plg[1][2] * p[12] * cd14 * swc[5];
        t[6] = (plg[1][1] * p[3] + plg[1][3] * p[4] + t71) * work.ctloc
            + (plg[1][1] * p[6] + plg[1][3] * p[7] + t72) * work.stloc;
    }

    // Semidiurnal.
    if sw[8] != 0.0 {
        let t81 = (plg[2][3] * p[23] + plg[2][5] * p[35]) * cd14 * swc[5];
        let t82 = (plg[2][3] * p[33] + plg[2][5] * p[36]) * cd14 * swc[5];
        t[7] = (plg[2][2] * p[5] + plg[2][4] * p[41] + t81) * work.c2tloc
            + (plg[2][2] * p[8] + plg[2][4] * p[42] + t82) * work.s2tloc;
    }

    // Terdiurnal.
    if sw[14] != 0.0 {
        t[13] = plg[3][3] * p[39] * work.s3tloc + plg[3][3] * p[40] * work.c3tloc;
    }

    // Magnetic activity.
    let sw9 = input.magnetic.sw();
    if sw9 != 0.0 {
        if sw9 == 1.0 {
            t[8] = work.apdf * (plg[0][2] * p[45] * swc[2] + p[32]);
        }
        if sw9 == -1.0 {
            t[8] = work.apt[0] * p[50] + plg[0][2] * p[96] * work.apt[0] * swc[2];
        }
    }

    // Longitudinal.
    if !(sw[10] == 0.0 || sw[11] == 0.0 || input.longitude_deg.re() <= -1000.0) {
        t[10] = (T::lift(1.0)
            + plg[0][1]
                * (p[80] * swc[5] * shifted(work.doy1, phase[P81], DR, doy, p[81])
                    + p[85] * swc[6] * shifted(work.doy2, phase[P86], 2.0 * DR, doy, p[86]))
            + p[83] * swc[3] * shifted(work.doy1, phase[P84], DR, doy, p[84])
            + p[87] * swc[4] * shifted(work.doy2, phase[P88], 2.0 * DR, doy, p[88]))
            * ((plg[1][2] * p[64]
                + plg[1][4] * p[65]
                + plg[1][6] * p[66]
                + plg[1][1] * p[74]
                + plg[1][3] * p[75]
                + plg[1][5] * p[76])
                * shifted(work.east, UNSHIFTED, DGTR, input.longitude_deg, 0.0)
                + (plg[1][2] * p[90]
                    + plg[1][4] * p[91]
                    + plg[1][6] * p[92]
                    + plg[1][1] * p[77]
                    + plg[1][3] * p[78]
                    + plg[1][5] * p[79])
                    * shifted_sin(work.east, UNSHIFTED, DGTR, input.longitude_deg, 0.0));
    }

    let mut tt = T::lift(0.0);
    for (index, term) in t.iter().enumerate() {
        let switch = if index == 8 { sw9 } else { sw[index + 1] };
        tt += *term * switch.abs();
    }
    tt
}

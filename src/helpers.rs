//! The small pieces the profiles are built from.
//!
//! Each is a direct translation of the like-named function in the reference, and each is
//! individually checkable, which the large ones are not.
//!
//! # Why the expressions look the way they do
//!
//! Floating-point addition and multiplication are not associative, so matching the
//! reference to the last bit means performing its operations in its order. Two habits
//! follow from that, and both are deliberate wherever they look awkward:
//!
//! - A constant that the reference writes on the *left* of an operator stays there, as
//!   `T::lift(3.0) / dx` rather than `dx.recip() * 3.0`. Commuting a product or reordering
//!   a sum is exact and is done freely; re-associating one is not and is not done.
//! - `pow(x, 2.0)` is translated as `powf(2.0)` and not as `x * x`. The two agree to within
//!   a rounding of the last bit and almost certainly agree exactly, but "almost certainly"
//!   is not what the fixture checks, and this is the faithful pass. Replacing it is a
//!   measured change for later.

use crate::scalar::Scalar;
use crate::work::Work;

/// The gas constant the model uses.
pub(crate) const RGAS: f64 = 831.4;

/// Degrees to radians, to the digits the reference carries.
///
/// It is not `PI / 180.0`, which differs in the sixth digit. Using a rounder value moves
/// every answer the model gives, so this is the reference's own constant.
pub(crate) const DGTR: f64 = 1.74533e-2;

/// Radians per day of year, for the annual terms.
pub(crate) const DR: f64 = 1.72142e-2;

/// Radians per second of universal time, for the UT terms.
pub(crate) const SR: f64 = 7.2722e-5;

/// Radians per hour of local solar time, for the diurnal terms.
pub(crate) const HR: f64 = 0.2618;

/// The most spline nodes any caller asks for.
pub(crate) const MAX_NODES: usize = 5;

/// Gravity and the effective Earth radius at a latitude.
///
/// `glatf`: gravity in cm/s^2, and a radius in km chosen so that the geopotential height
/// of [`zeta`] is right.
pub(crate) fn glatf<T: Scalar>(lat: T) -> (T, T) {
    let c2 = (lat * (2.0 * DGTR)).cos();
    let gv = (c2 * -0.0026373 + 1.0) * 980.616;
    let reff = (gv * 2.0) / (c2 * 2.27e-9 + 3.085462e-6) * 1.0e-5;
    (gv, reff)
}

/// Geopotential height of `zz` above `zl`, in km.
///
/// `zeta`. The model's vertical coordinate: an altitude difference scaled so that a
/// hydrostatic profile computed with constant gravity is right in a field that falls off.
#[inline]
pub(crate) fn zeta<T: Scalar>(work: &Work<T>, zz: T, zl: T) -> T {
    (zz - zl) * (work.re + zl) / (work.re + zz)
}

/// Chemistry and dissociation correction.
///
/// `ccor`, a multiplicative factor on a species' density. `r` is the target ratio, `h1` the
/// transition scale length and `zh` the altitude at which half of `r` has been applied.
pub(crate) fn ccor<T: Scalar>(alt: T, r: T, h1: f64, zh: f64) -> T {
    let e = (alt - zh) / h1;
    // The reference returns exp(0) and exp(r) at the ends rather than clamping the
    // exponent, so the factor goes to 1 above the transition and to exp(r) below it. Both
    // ends are constant in the altitude, and so have no derivative.
    if e.re() > 70.0 {
        return T::lift(1.0);
    }
    if e.re() < -70.0 {
        return r.exp();
    }
    (r / (e.exp() + 1.0)).exp()
}

/// Chemistry and dissociation correction with two transition scale lengths.
///
/// `ccor2`. As [`ccor`], but the transition is the mean of two of different widths.
pub(crate) fn ccor2<T: Scalar>(alt: T, r: T, h1: f64, zh: f64, h2: f64) -> T {
    let e1 = (alt - zh) / h1;
    let e2 = (alt - zh) / h2;
    if e1.re() > 70.0 || e2.re() > 70.0 {
        return T::lift(1.0);
    }
    if e1.re() < -70.0 && e2.re() < -70.0 {
        return r.exp();
    }
    (r / ((e1.exp() + e2.exp()) * 0.5 + 1.0)).exp()
}

/// The scale height of a species of mass `xm` at a temperature, in km.
///
/// `scalh`.
pub(crate) fn scalh<T: Scalar>(work: &Work<T>, alt: f64, xm: f64, temp: T) -> T {
    let base = T::lift(alt) / work.re + 1.0;
    let g = work.gsurf / (base * base);
    temp * RGAS / (g * xm)
}

/// The turbopause correction: a root-mean blend of a diffusive and a mixed density.
///
/// `dnet`. `zhm` is the transition scale length, `xmm` the fully mixed molecular weight and
/// `xm` the species weight.
pub(crate) fn dnet<T: Scalar>(dd: T, dm: T, zhm: f64, xmm: f64, xm: f64) -> T {
    let a = zhm / (xmm - xm);
    if !(dm.re() > 0.0 && dd.re() > 0.0) {
        // The reference prints a message here and carries on with whichever of the two it
        // can use. Reaching it means a density came out non-positive, which the blend has
        // no logarithm for; it happens when a switch has turned off a term the blend needs.
        if dd.re() == 0.0 && dm.re() == 0.0 {
            return T::lift(1.0);
        }
        if dm.re() == 0.0 {
            return dd;
        }
        if dd.re() == 0.0 {
            return dm;
        }
    }
    let ylog = (dm / dd).ln() * a;
    if ylog.re() < -10.0 {
        return dd;
    }
    if ylog.re() > 10.0 {
        return dm;
    }
    dd * (ylog.exp() + 1.0).powf(1.0 / a)
}

/// The second derivatives of a cubic spline through `(x, y)`.
///
/// `spline`, which is Numerical Recipes' routine of that name. `yp1` and `ypn` are the
/// first derivatives imposed at the ends; a value at or above 1e30 would ask for a zero
/// second derivative there instead, which this model never does.
pub(crate) fn spline<T: Scalar>(x: &[T], y: &[T], yp1: T, ypn: T, y2: &mut [T]) {
    let n = x.len();
    // The reference heap-allocates this per call. Five nodes is the most any caller wants.
    let mut u = [T::lift(0.0); MAX_NODES];

    if yp1.re() > 0.99e30 {
        y2[0] = T::lift(0.0);
        u[0] = T::lift(0.0);
    } else {
        y2[0] = T::lift(-0.5);
        u[0] = (T::lift(3.0) / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for i in 1..n - 1 {
        let sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        let p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
            * 6.0
            / (x[i + 1] - x[i - 1])
            - sig * u[i - 1])
            / p;
    }
    let (qn, un) = if ypn.re() > 0.99e30 {
        (T::lift(0.0), T::lift(0.0))
    } else {
        (
            T::lift(0.5),
            (T::lift(3.0) / (x[n - 1] - x[n - 2]))
                * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])),
        )
    };
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for k in (0..n - 1).rev() {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
}

/// A cubic spline evaluated at `x`.
///
/// `splint`. `y2a` is what [`spline`] produced for the same nodes.
pub(crate) fn splint<T: Scalar>(xa: &[T], ya: &[T], y2a: &[T], x: T) -> T {
    let mut klo = 0usize;
    let mut khi = xa.len() - 1;
    while khi - klo > 1 {
        let k = (khi + klo) / 2;
        if xa[k].re() > x.re() {
            khi = k;
        } else {
            klo = k;
        }
    }
    let h = xa[khi] - xa[klo];
    let a = (xa[khi] - x) / h;
    let b = (x - xa[klo]) / h;
    a * ya[klo]
        + b * ya[khi]
        + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0
}

/// A cubic spline integrated from `xa[0]` to `x`.
///
/// `splini`. The integral of the interpolant, which is what turns a profile of inverse
/// temperature into a density ratio.
pub(crate) fn splini<T: Scalar>(xa: &[T], ya: &[T], y2a: &[T], x: T) -> T {
    let n = xa.len();
    let mut yi = T::lift(0.0);
    let mut klo = 0usize;
    let mut khi = 1usize;
    while x.re() > xa[klo].re() && khi < n {
        // Every interval below the last is integrated whole; the one holding x is
        // integrated up to x.
        //
        // The reference's test is `x < xa[khi]`, which takes the node when the two are
        // equal. That is the same number, so the value is unchanged, but at a dual type it
        // is the difference between carrying dx and dropping it: substituting the node
        // replaces x by a constant and the integral comes out with no derivative at all.
        // Sampling exactly on a node is not a corner case here — the mesosphere's nodes
        // are at 0, 10, 15, 20 and 32.5 km, which are round numbers a caller will ask for.
        let mut xx = x;
        if khi < n - 1 && x.re() > xa[khi].re() {
            xx = xa[khi];
        }
        let h = xa[khi] - xa[klo];
        let a = (xa[khi] - xx) / h;
        let b = (xx - xa[klo]) / h;
        let a2 = a * a;
        let b2 = b * b;
        yi += ((T::lift(1.0) - a2) * ya[klo] / 2.0
            + b2 * ya[khi] / 2.0
            + ((-(a2 * a2 + 1.0) / 4.0 + a2 / 2.0) * y2a[klo]
                + (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi])
                * h
                * h
                / 6.0)
            * h;
        klo += 1;
        khi += 1;
    }
    yi
}

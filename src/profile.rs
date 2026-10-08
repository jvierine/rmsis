//! The temperature and density profiles: `densu` above 72.5 km and `densm` below it.
//!
//! Both work the same way. A spline is laid through the inverse of the temperature at a
//! handful of nodes, because that is what makes the hydrostatic integral a spline integral
//! and so closed-form; the temperature at the altitude asked for comes from evaluating it,
//! and the density from integrating it.
//!
//! The reference passes the node arrays in as pointers and every caller passes the same
//! ones, so here they are read from [`Work`] directly: `densu` always gets `tn1` and
//! `tgn1`, `densm` always gets `tn2`/`tgn2` and then `tn3`/`tgn3`.

use crate::helpers::{spline, splini, splint, zeta, MAX_NODES, RGAS};
use crate::scalar::Scalar;
use crate::work::Work;

/// Temperature and density above the lower thermosphere's first node.
///
/// `densu`. Returns the density of a species of mass `xm`, or the temperature when `xm` is
/// zero, and writes the temperature at the altitude to `tz` either way.
///
/// `dlb` and `tlb` are the density and temperature at the lower boundary `zlb`, `tinf` the
/// exospheric temperature, `alpha` the species' thermal diffusion coefficient and `s2` the
/// temperature gradient parameter.
///
/// It also writes the Bates profile's temperature and gradient into `tn1[0]` and `tgn1[0]`,
/// which is how the spline below `za` is joined continuously to the profile above it. The
/// mesosphere reads those afterwards, so the order of calls matters.
#[allow(clippy::too_many_arguments)]
pub(crate) fn densu<T: Scalar>(
    work: &mut Work<T>,
    alt: T,
    dlb: T,
    tinf: T,
    tlb: T,
    xm: f64,
    alpha: f64,
    tz: &mut T,
    zlb: f64,
    s2: T,
    zn1: &[f64; MAX_NODES],
) -> T {
    let za = zn1[0];
    let z = if alt.re() > za { alt } else { T::lift(za) };

    // Geopotential height above the lower boundary, and the Bates temperature there.
    let zg2 = zeta(work, z, T::lift(zlb));
    let tt = tinf - (tinf - tlb) * ((-s2) * zg2).exp();
    let ta = tt;
    *tz = tt;
    let mut densu_temp = *tz;

    // Carried out of the branch below, because the density part reads them.
    let mut xs = [T::lift(0.0); MAX_NODES];
    let mut ys = [T::lift(0.0); MAX_NODES];
    let mut y2out = [T::lift(0.0); MAX_NODES];
    let mut x = T::lift(0.0);
    let mut zgdif = T::lift(0.0);
    let mut z1 = 0.0f64;
    let mut t1 = T::lift(0.0);

    if alt.re() < za {
        // The temperature gradient at za from the Bates profile, which the spline below it
        // is made to match.
        let ratio = (work.re + zlb) / (work.re + za);
        let dta = (tinf - ta) * s2 * (ratio * ratio);
        work.tgn1[0] = dta;
        work.tn1[0] = ta;

        let mn = zn1.len();
        let z = if alt.re() > zn1[mn - 1] {
            alt
        } else {
            T::lift(zn1[mn - 1])
        };
        z1 = zn1[0];
        let z2 = zn1[mn - 1];
        t1 = work.tn1[0];
        let t2 = work.tn1[mn - 1];
        let zg = zeta(work, z, T::lift(z1));
        zgdif = zeta(work, T::lift(z2), T::lift(z1));

        for k in 0..mn {
            xs[k] = zeta(work, T::lift(zn1[k]), T::lift(z1)) / zgdif;
            ys[k] = T::lift(1.0) / work.tn1[k];
        }
        let yd1 = -work.tgn1[0] / (t1 * t1) * zgdif;
        let yd2 = -work.tgn1[1] / (t2 * t2) * zgdif * ((work.re + z2) / (work.re + z1)).powf(2.0);

        spline(&xs[..mn], &ys[..mn], yd1, yd2, &mut y2out[..mn]);
        x = zg / zgdif;
        *tz = T::lift(1.0) / splint(&xs[..mn], &ys[..mn], &y2out[..mn], x);
        densu_temp = *tz;
    }
    if xm == 0.0 {
        return densu_temp;
    }

    // Density above za, from the Bates profile.
    let lower = T::lift(zlb) / work.re + 1.0;
    let glb = work.gsurf / (lower * lower);
    let gamma = glb * xm / (s2 * RGAS * tinf);
    let mut expl = ((-s2) * gamma * zg2).exp();
    if expl.re() > 50.0 {
        expl = T::lift(50.0);
    }
    if tt.re() <= 0.0 {
        expl = T::lift(50.0);
    }
    densu_temp = dlb * (tlb / tt).powd(gamma + (1.0 + alpha)) * expl;
    if alt.re() >= za {
        return densu_temp;
    }

    // Density below za, from the integral of the spline.
    let below = T::lift(z1) / work.re + 1.0;
    let glb = work.gsurf / (below * below);
    let gamm = glb * xm * zgdif / RGAS;
    let mn = zn1.len();
    let yi = splini(&xs[..mn], &ys[..mn], &y2out[..mn], x);
    let mut expl = gamm * yi;
    if expl.re() > 50.0 {
        expl = T::lift(50.0);
    }
    if tz.re() <= 0.0 {
        expl = T::lift(50.0);
    }
    densu_temp * (t1 / *tz).powf(1.0 + alpha) * (-expl).exp()
}

/// Temperature and density below 72.5 km.
///
/// `densm`. Returns the density of a fully mixed atmosphere of mass `xm` scaled from `d0`,
/// or the temperature when `xm` is zero, and writes the temperature to `tz` either way.
///
/// Two layers, stratosphere and mesosphere down to 32.5 km and then troposphere below it,
/// each a spline through its own nodes.
pub(crate) fn densm<T: Scalar>(
    work: &Work<T>,
    alt: T,
    d0: T,
    xm: f64,
    tz: &mut T,
    zn3: &[f64; 5],
    zn2: &[f64; 4],
) -> T {
    let mut densm_tmp = d0;

    if alt.re() > zn2[0] {
        // Above where this routine applies. The reference returns the temperature it was
        // handed, which it has not written; no caller reaches here.
        return if xm == 0.0 { *tz } else { d0 };
    }

    // Stratosphere and mesosphere.
    let z = if alt.re() > zn2[zn2.len() - 1] {
        alt
    } else {
        T::lift(zn2[zn2.len() - 1])
    };
    densm_tmp = layer(work, z, zn2, &work.tn2, &work.tgn2, xm, densm_tmp, tz);

    if alt.re() > zn3[0] {
        return if xm == 0.0 { *tz } else { densm_tmp };
    }

    // Troposphere and lower stratosphere. The altitude is used as it is here, not clamped.
    let densm_tmp = layer(work, alt, zn3, &work.tn3, &work.tgn3, xm, densm_tmp, tz);
    if xm == 0.0 {
        *tz
    } else {
        densm_tmp
    }
}

/// One spline layer of [`densm`]: the temperature at `z`, and the density scaled across it.
///
/// The reference writes this out twice, once per layer, with the node arrays swapped.
#[allow(clippy::too_many_arguments)]
fn layer<T: Scalar>(
    work: &Work<T>,
    z: T,
    zn: &[f64],
    tn: &[T],
    tgn: &[T; 2],
    xm: f64,
    density: T,
    tz: &mut T,
) -> T {
    let mn = zn.len();
    let z1 = zn[0];
    let z2 = zn[mn - 1];
    let t1 = tn[0];
    let t2 = tn[mn - 1];
    let zg = zeta(work, z, T::lift(z1));
    let zgdif = zeta(work, T::lift(z2), T::lift(z1));

    let mut xs = [T::lift(0.0); MAX_NODES];
    let mut ys = [T::lift(0.0); MAX_NODES];
    let mut y2out = [T::lift(0.0); MAX_NODES];
    for k in 0..mn {
        xs[k] = zeta(work, T::lift(zn[k]), T::lift(z1)) / zgdif;
        ys[k] = T::lift(1.0) / tn[k];
    }
    let yd1 = -tgn[0] / (t1 * t1) * zgdif;
    let edge = (work.re + z2) / (work.re + z1);
    let yd2 = -tgn[1] / (t2 * t2) * zgdif * (edge * edge);

    spline(&xs[..mn], &ys[..mn], yd1, yd2, &mut y2out[..mn]);
    let x = zg / zgdif;
    *tz = T::lift(1.0) / splint(&xs[..mn], &ys[..mn], &y2out[..mn], x);

    if xm == 0.0 {
        return density;
    }

    let below = T::lift(z1) / work.re + 1.0;
    let glb = work.gsurf / (below * below);
    let gamm = glb * xm * zgdif / RGAS;
    let yi = splini(&xs[..mn], &ys[..mn], &y2out[..mn], x);
    let mut expl = gamm * yi;
    if expl.re() > 50.0 {
        expl = T::lift(50.0);
    }
    density * (t1 / *tz) * (-expl).exp()
}

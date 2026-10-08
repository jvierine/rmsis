//! The model itself: `gts7` for the thermosphere, and `gtd7` and `gtd7d` over the whole
//! range.

use crate::expansion::{glob7s, globe7, prepare, Phase, PHASES};
use crate::helpers::{ccor, ccor2, dnet, glatf, scalh, DGTR, DR, MAX_NODES};
use crate::input::{Input, Output};
use crate::profile::{densm, densu};
use crate::scalar::Scalar;
use crate::switches::{Resolved, Switches, Units};
use crate::tables::{PAVGM, PD, PDL, PDM, PMA, PS, PT, PTL, PTM};
use crate::work::Work;

/// Thermal diffusion coefficients, by the species order of the output.
///
/// `alpha` in the reference. Only helium, argon and hydrogen have one.
const ALPHA: [f64; 9] = [-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0];

/// Altitudes above which each species' mixed correction is dropped, km.
///
/// `altl` in the reference, in the order N2, He, O, O2, Ar, H, N — the reference's own
/// order here, which is not the output's.
const ALTL: [f64; 8] = [200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0];

/// Nodes of the troposphere and lower stratosphere temperature spline, km.
const ZN3: [f64; 5] = [32.5, 20.0, 15.0, 10.0, 0.0];

/// Nodes of the upper stratosphere and mesosphere temperature spline, km.
const ZN2: [f64; 4] = [72.5, 55.0, 45.0, 32.5];

/// Altitude above which mixing is no longer complete, km. `zmix` in the reference.
const ZMIX: f64 = 62.5;

/// Mass of a nucleon in g, as the reference writes it, for summing a mass density.
const NUCLEON_MASS: f64 = 1.66e-24;

/// NRLMSISE-00.
///
/// Built once from a set of [`Switches`] and then evaluated as often as wanted. Building it
/// resolves the switches into the form the expansions multiply by, which is the reference's
/// `tselec` done once rather than on every call.
///
/// ```
/// use rmsis::{Input, MagneticActivity, Model};
///
/// let model = Model::new();
/// let output = model.gtd7d::<f64>(&Input::consistent(
///     172,                                     // day of year
///     29000.0,                                 // seconds of universal time
///     400.0,                                   // km
///     60.0,                                    // degrees north
///     -70.0,                                   // degrees east
///     150.0,                                   // 81-day average F10.7
///     150.0,                                   // previous day's F10.7
///     MagneticActivity::Daily { ap: 4.0 },
/// ));
/// // A couple of picograms per cubic meter at 400 km, in SI by default.
/// assert!((output.total_mass - 2.4e-12).abs() < 0.1e-12);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    switches: Switches,
    resolved: Resolved,
    phases: Phases,
}

/// The phase shifts of every coefficient table, worked out once when a model is built.
///
/// Each entry is the sine and cosine of one table's day-of-year phase coefficients, which
/// is what lets the `fast` path reach a phase-shifted cosine by an angle-sum rather than by
/// calling `cos`. Without that feature a [`Phase`] is zero-sized and this costs nothing.
#[derive(Clone, Debug)]
struct Phases {
    pt: [Phase; PHASES],
    ps: [Phase; PHASES],
    pd: [[Phase; PHASES]; 9],
    ptl: [[Phase; PHASES]; 4],
    pma: [[Phase; PHASES]; 10],
}

impl Phases {
    fn new() -> Self {
        Self {
            pt: Phase::of(&PT),
            ps: Phase::of(&PS),
            pd: core::array::from_fn(|row| Phase::of(&PD[row])),
            ptl: core::array::from_fn(|row| Phase::of(&PTL[row])),
            pma: core::array::from_fn(|row| Phase::of(&PMA[row])),
        }
    }
}

impl Default for Model {
    fn default() -> Self {
        Self::new()
    }
}

impl Model {
    /// The model with every variation included, answering in SI.
    pub fn new() -> Self {
        Self::with_switches(Switches::default())
    }

    /// The model with a chosen set of variations and units.
    pub fn with_switches(switches: Switches) -> Self {
        Self {
            switches,
            resolved: switches.resolve(),
            phases: Phases::new(),
        }
    }

    /// Which variations this model includes.
    pub fn switches(&self) -> Switches {
        self.switches
    }

    /// Density and temperature from the surface to the lower exosphere.
    ///
    /// The total mass density sums helium, atomic oxygen, molecular nitrogen, molecular
    /// oxygen, argon, hydrogen and atomic nitrogen. It leaves out anomalous oxygen, which
    /// is what [`Model::gtd7d`] adds, and a drag calculation wants that one. The individual
    /// number densities are the same either way, anomalous oxygen included.
    pub fn gtd7<T: Scalar>(&self, input: &Input<T>) -> Output<T> {
        let mut work = Work::new();
        self.evaluate(input, &mut work)
    }

    /// Density and temperature, with anomalous oxygen in the total mass density.
    ///
    /// This is the total a drag calculation wants: the reference calls it the effective
    /// total mass density for drag. Anomalous oxygen is hot atomic oxygen, and above about
    /// 500 km it is a large part of the mass that a satellite runs into.
    ///
    /// Nothing but [`Output::total_mass`] differs from [`Model::gtd7`].
    pub fn gtd7d<T: Scalar>(&self, input: &Input<T>) -> Output<T> {
        let mut work = Work::new();
        let mut output = self.evaluate(input, &mut work);
        output.total_mass = mass_density(&output, true);
        if self.switches.units == Units::Si {
            output.total_mass /= 1000.0;
        }
        output
    }

    /// Density and temperature from the thermospheric part of the model alone.
    ///
    /// `gts7`. Only meaningful above 72.5 km; [`Model::gtd7`] handles the whole range, and
    /// below 72.5 km it calls this at 72.5 and blends downward from there.
    pub fn gts7<T: Scalar>(&self, input: &Input<T>) -> Output<T> {
        let mut work = Work::new();
        let xlat = if self.resolved.sw[2] == 0.0 {
            T::lift(45.0)
        } else {
            input.latitude_deg
        };
        (work.gsurf, work.re) = glatf(xlat);
        prepare(input, &self.resolved, &mut work);
        self.thermosphere(input, input.altitude_km, &mut work)
    }

    /// `gtd7`: the thermosphere, and below 72.5 km the mesosphere blended onto it.
    fn evaluate<T: Scalar>(&self, input: &Input<T>, work: &mut Work<T>) -> Output<T> {
        let (sw, metric) = (&self.resolved.sw, self.switches.units == Units::Si);

        // Latitude variation of gravity, which switch 2 also governs.
        let xlat = if sw[2] == 0.0 {
            T::lift(45.0)
        } else {
            input.latitude_deg
        };
        (work.gsurf, work.re) = glatf(xlat);
        // The parts of the expansions that depend only on the input, computed once for
        // the whole evaluation rather than on each of globe7's eleven calls.
        prepare(input, &self.resolved, work);

        let xmm = PDM[2][4];

        // The thermosphere, evaluated no lower than the top of the mesosphere spline.
        let altt = if input.altitude_km.re() > ZN2[0] {
            input.altitude_km
        } else {
            T::lift(ZN2[0])
        };
        let thermosphere = self.thermosphere(input, altt, work);

        let dm28m = if metric { work.dm28 * 1.0e6 } else { work.dm28 };
        let mut output = thermosphere;

        if input.altitude_km.re() >= ZN2[0] {
            return output;
        }

        // Lower mesosphere and upper stratosphere, 32.5 to 72.5 km. The inverse of the
        // temperature at each node is a linear function of the spherical harmonics.
        work.tgn2[0] = work.tgn1[1];
        work.tn2[0] = work.tn1[4];
        work.tn2[1] = T::lift(PMA[0][0] * PAVGM[0])
            / (-(glob7s(&PMA[0], &self.phases.pma[0], input, &self.resolved, work) * sw[20]) + 1.0);
        work.tn2[2] = T::lift(PMA[1][0] * PAVGM[1])
            / (-(glob7s(&PMA[1], &self.phases.pma[1], input, &self.resolved, work) * sw[20]) + 1.0);
        work.tn2[3] = T::lift(PMA[2][0] * PAVGM[2])
            / (-(glob7s(&PMA[2], &self.phases.pma[2], input, &self.resolved, work)
                * sw[20]
                * sw[22])
                + 1.0);
        let tgn =
            (glob7s(&PMA[9], &self.phases.pma[9], input, &self.resolved, work) * sw[20] * sw[22]
                + 1.0)
                * (PAVGM[8] * PMA[9][0])
                * work.tn2[3]
                * work.tn2[3]
                / (PMA[2][0] * PAVGM[2]).powi(2);
        work.tgn2[1] = tgn;
        work.tn3[0] = work.tn2[3];

        // Lower stratosphere and troposphere, at and below 32.5 km.
        //
        // The reference's test here is `alt < 32.5`, while `densm` enters this layer at
        // `alt <= 32.5`. At exactly 32.5 km it therefore splines through nodes its own
        // call never filled, and reads whichever values the previous call to the model
        // left in them. In C those are file-scope statics, so the answer there depends on
        // what was asked before; here they would be zero, and the spline's second
        // derivatives would come out infinite, and the density NaN.
        //
        // These nodes do not depend on the altitude — they are the expansions in day,
        // position and flux — so filling them at exactly 32.5 km gives precisely the
        // values a preceding lower call would have left. That makes this agree with the
        // reference's numbers wherever those numbers mean anything, which the grid fixture
        // checks at 32.5 km in sixty scenarios, and removes the NaN. It is the one place
        // the port departs from the reference's control flow.
        if input.altitude_km.re() <= ZN3[0] {
            work.tgn3[0] = work.tgn2[1];
            work.tn3[1] = T::lift(PMA[3][0] * PAVGM[3])
                / (-(glob7s(&PMA[3], &self.phases.pma[3], input, &self.resolved, work) * sw[22])
                    + 1.0);
            work.tn3[2] = T::lift(PMA[4][0] * PAVGM[4])
                / (-(glob7s(&PMA[4], &self.phases.pma[4], input, &self.resolved, work) * sw[22])
                    + 1.0);
            work.tn3[3] = T::lift(PMA[5][0] * PAVGM[5])
                / (-(glob7s(&PMA[5], &self.phases.pma[5], input, &self.resolved, work) * sw[22])
                    + 1.0);
            work.tn3[4] = T::lift(PMA[6][0] * PAVGM[6])
                / (-(glob7s(&PMA[6], &self.phases.pma[6], input, &self.resolved, work) * sw[22])
                    + 1.0);
            let tgn = (glob7s(&PMA[7], &self.phases.pma[7], input, &self.resolved, work) * sw[22]
                + 1.0)
                * (PMA[7][0] * PAVGM[7])
                * work.tn3[4]
                * work.tn3[4]
                / (PMA[6][0] * PAVGM[6]).powi(2);
            work.tgn3[1] = tgn;
        }

        // Linear transition to full mixing below 72.5 km.
        let mut dmc = T::lift(0.0);
        if input.altitude_km.re() > ZMIX {
            dmc = -((T::lift(ZN2[0]) - input.altitude_km) / (ZN2[0] - ZMIX)) + 1.0;
        }
        let dz28 = thermosphere.molecular_nitrogen;
        let mut tz = T::lift(0.0);

        // Molecular nitrogen, which the others are scaled from.
        let dmr = thermosphere.molecular_nitrogen / dm28m - 1.0;
        let dm = densm(work, input.altitude_km, dm28m, xmm, &mut tz, &ZN3, &ZN2);
        output.molecular_nitrogen = dm * (dmr * dmc + 1.0);

        // Helium.
        let dmr = thermosphere.helium / (dz28 * PDM[0][1]) - 1.0;
        output.helium = output.molecular_nitrogen * PDM[0][1] * (dmr * dmc + 1.0);

        // Atomic oxygen and anomalous oxygen are not carried below 72.5 km.
        output.atomic_oxygen = T::lift(0.0);
        output.anomalous_oxygen = T::lift(0.0);

        // Molecular oxygen.
        let dmr = thermosphere.molecular_oxygen / (dz28 * PDM[3][1]) - 1.0;
        output.molecular_oxygen = output.molecular_nitrogen * PDM[3][1] * (dmr * dmc + 1.0);

        // Argon.
        let dmr = thermosphere.argon / (dz28 * PDM[4][1]) - 1.0;
        output.argon = output.molecular_nitrogen * PDM[4][1] * (dmr * dmc + 1.0);

        // Hydrogen and atomic nitrogen are not carried below 72.5 km either.
        output.hydrogen = T::lift(0.0);
        output.atomic_nitrogen = T::lift(0.0);

        // The total is summed from the densities as they now stand, which are already in
        // the output units, so only the thousand is applied here.
        output.total_mass = mass_density(&output, false);
        if metric {
            output.total_mass /= 1000.0;
        }

        // The temperature at the altitude, from the same profile with no mass.
        densm(
            work,
            input.altitude_km,
            T::lift(1.0),
            0.0,
            &mut tz,
            &ZN3,
            &ZN2,
        );
        output.temperature = tz;
        output
    }

    /// The thermospheric part: `gts7`.
    ///
    /// Each species is built the same way: a variation factor at the lower boundary, a
    /// diffusive density from [`densu`], and then, where switch 15 allows it and the
    /// altitude is low enough, a fully mixed density blended in by [`dnet`] and one or two
    /// chemistry corrections from [`ccor`] and [`ccor2`].
    fn thermosphere<T: Scalar>(&self, input: &Input<T>, alt: T, work: &mut Work<T>) -> Output<T> {
        let sw = &self.resolved.sw;
        let mut output = Output::zeroed();

        // The first node of the lower thermosphere spline is itself a fitted coefficient.
        let za = PDL[1][15];
        let zn1: [f64; MAX_NODES] = [za, 110.0, 100.0, 90.0, 72.5];

        // Variation in the exospheric temperature does not matter below za.
        let tinf = if alt.re() > zn1[0] {
            T::lift(PTM[0] * PT[0])
                * (globe7(&PT, &self.phases.pt, input, &self.resolved, work) * sw[16] + 1.0)
        } else {
            T::lift(PTM[0] * PT[0])
        };
        output.exospheric_temperature = tinf;

        // Variation in the 120 km gradient does not matter below the bottom node.
        let g0 = if alt.re() > zn1[4] {
            T::lift(PTM[3] * PS[0])
                * (globe7(&PS, &self.phases.ps, input, &self.resolved, work) * sw[19] + 1.0)
        } else {
            T::lift(PTM[3] * PS[0])
        };
        let tlb = (globe7(&PD[3], &self.phases.pd[3], input, &self.resolved, work) * sw[17] + 1.0)
            * PTM[1]
            * PD[3][0];
        let s = g0 / (tinf - tlb);

        // Lower thermosphere temperatures do not matter for density above 300 km.
        if alt.re() < 300.0 {
            work.tn1[1] = T::lift(PTM[6] * PTL[0][0])
                / (-(glob7s(&PTL[0], &self.phases.ptl[0], input, &self.resolved, work) * sw[18])
                    + 1.0);
            work.tn1[2] = T::lift(PTM[2] * PTL[1][0])
                / (-(glob7s(&PTL[1], &self.phases.ptl[1], input, &self.resolved, work) * sw[18])
                    + 1.0);
            work.tn1[3] = T::lift(PTM[7] * PTL[2][0])
                / (-(glob7s(&PTL[2], &self.phases.ptl[2], input, &self.resolved, work) * sw[18])
                    + 1.0);
            work.tn1[4] = T::lift(PTM[4] * PTL[3][0])
                / (-(glob7s(&PTL[3], &self.phases.ptl[3], input, &self.resolved, work)
                    * sw[18]
                    * sw[20])
                    + 1.0);
            let tgn = (glob7s(&PMA[8], &self.phases.pma[8], input, &self.resolved, work)
                * sw[18]
                * sw[20]
                + 1.0)
                * (PTM[8] * PMA[8][0])
                * work.tn1[4]
                * work.tn1[4]
                / (PTM[4] * PTL[3][0]).powi(2);
            work.tgn1[1] = tgn;
        } else {
            work.tn1[1] = T::lift(PTM[6] * PTL[0][0]);
            work.tn1[2] = T::lift(PTM[2] * PTL[1][0]);
            work.tn1[3] = T::lift(PTM[7] * PTL[2][0]);
            work.tn1[4] = T::lift(PTM[4] * PTL[3][0]);
            work.tgn1[1] = T::lift(PTM[8] * PMA[8][0]) * work.tn1[4] * work.tn1[4]
                / (PTM[4] * PTL[3][0]).powi(2);
        }

        // Variation of the molecular nitrogen density at the lower boundary.
        let g28 = globe7(&PD[2], &self.phases.pd[2], input, &self.resolved, work) * sw[21];

        // Variation of the turbopause height.
        let zhf = ((input.latitude_deg * DGTR).sin()
            * (sw[5] * PDL[0][24])
            * (DR * (f64::from(input.day_of_year) - PT[13])).cos()
            + 1.0)
            * PDL[1][24];
        output.exospheric_temperature = tinf;
        let xmm = PDM[2][4];
        let z = alt;
        let mut tz = T::lift(0.0);

        // Molecular nitrogen, which every other species below 72.5 km is scaled from.
        let db28 = g28.exp() * PDM[2][0] * PD[2][0];
        output.molecular_nitrogen = densu(
            work, z, db28, tinf, tlb, 28.0, ALPHA[2], &mut tz, PTM[5], s, &zn1,
        );
        let zh28 = zhf * PDM[2][2];
        let zhm28 = PDM[2][3] * PDL[1][5];
        let xmd = 28.0 - xmm;
        let b28 = densu(
            work,
            zh28,
            db28,
            tinf,
            tlb,
            xmd,
            ALPHA[2] - 1.0,
            &mut tz,
            PTM[5],
            s,
            &zn1,
        );
        // The reference computes the mixed density only when switch 15 is on, but `gtd7`
        // reads it below 72.5 km whatever the switch says, where it is the density every
        // other species is scaled from. With the switch off the reference therefore reads
        // a value from whichever call last wrote it — a file-scope static — and returns an
        // answer that depends on what was asked before. Computing it regardless costs one
        // `densu` call in a configuration nobody runs by default, and changes nothing when
        // the switch is on: the blend below is still what the switch governs.
        if z.re() <= ALTL[2] {
            let dm28 = densu(
                work, z, b28, tinf, tlb, xmm, ALPHA[2], &mut tz, PTM[5], s, &zn1,
            );
            work.dm28 = dm28;
            if sw[15] != 0.0 {
                output.molecular_nitrogen = dnet(output.molecular_nitrogen, dm28, zhm28, xmm, 28.0);
            }
        }

        // Helium.
        let g4 = globe7(&PD[0], &self.phases.pd[0], input, &self.resolved, work) * sw[21];
        let db04 = g4.exp() * PDM[0][0] * PD[0][0];
        output.helium = densu(
            work, z, db04, tinf, tlb, 4.0, ALPHA[0], &mut tz, PTM[5], s, &zn1,
        );
        // Original Fortran uses Z.GT.ALTL to skip: include the boundary.
        if sw[15] != 0.0 && z.re() <= ALTL[0] {
            let zh04 = PDM[0][2];
            let b04 = densu(
                work,
                T::lift(zh04),
                db04,
                tinf,
                tlb,
                4.0 - xmm,
                ALPHA[0] - 1.0,
                &mut tz,
                PTM[5],
                s,
                &zn1,
            );
            let dm04 = densu(work, z, b04, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
            work.dm04 = dm04;
            let zhm04 = zhm28;
            output.helium = dnet(output.helium, dm04, zhm04, xmm, 4.0);
            // Correction to the specified mixing ratio at the ground.
            let rl = (b28 * PDM[0][1] / b04).ln();
            let zc04 = PDM[0][4] * PDL[1][0];
            let hc04 = PDM[0][5] * PDL[1][1];
            output.helium *= ccor(z, rl, hc04, zc04);
        }

        // Atomic oxygen.
        let g16 = globe7(&PD[1], &self.phases.pd[1], input, &self.resolved, work) * sw[21];
        let db16 = g16.exp() * PDM[1][0] * PD[1][0];
        output.atomic_oxygen = densu(
            work, z, db16, tinf, tlb, 16.0, ALPHA[1], &mut tz, PTM[5], s, &zn1,
        );
        if sw[15] != 0.0 && z.re() <= ALTL[1] {
            let zh16 = PDM[1][2];
            let b16 = densu(
                work,
                T::lift(zh16),
                db16,
                tinf,
                tlb,
                16.0 - xmm,
                ALPHA[1] - 1.0,
                &mut tz,
                PTM[5],
                s,
                &zn1,
            );
            let dm16 = densu(work, z, b16, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
            work.dm16 = dm16;
            let zhm16 = zhm28;
            output.atomic_oxygen = dnet(output.atomic_oxygen, dm16, zhm16, xmm, 16.0);
            let rl = ((input.f107_average - 150.0) * (sw[1] * PDL[0][23]) + 1.0)
                * (PDM[1][1] * PDL[1][16]);
            let hc16 = PDM[1][5] * PDL[1][3];
            let zc16 = PDM[1][4] * PDL[1][2];
            let hc216 = PDM[1][5] * PDL[1][4];
            output.atomic_oxygen *= ccor2(z, rl, hc16, zc16, hc216);
            // Chemistry correction.
            let hcc16 = PDM[1][7] * PDL[1][13];
            let zcc16 = PDM[1][6] * PDL[1][12];
            let rc16 = PDM[1][3] * PDL[1][14];
            output.atomic_oxygen *= ccor(z, T::lift(rc16), hcc16, zcc16);
        }

        // Molecular oxygen.
        let g32 = globe7(&PD[4], &self.phases.pd[4], input, &self.resolved, work) * sw[21];
        let db32 = g32.exp() * PDM[3][0] * PD[4][0];
        output.molecular_oxygen = densu(
            work, z, db32, tinf, tlb, 32.0, ALPHA[3], &mut tz, PTM[5], s, &zn1,
        );
        if sw[15] != 0.0 {
            if z.re() <= ALTL[3] {
                let zh32 = PDM[3][2];
                let b32 = densu(
                    work,
                    T::lift(zh32),
                    db32,
                    tinf,
                    tlb,
                    32.0 - xmm,
                    ALPHA[3] - 1.0,
                    &mut tz,
                    PTM[5],
                    s,
                    &zn1,
                );
                let dm32 = densu(work, z, b32, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
                work.dm32 = dm32;
                let zhm32 = zhm28;
                output.molecular_oxygen = dnet(output.molecular_oxygen, dm32, zhm32, xmm, 32.0);
                let rl = (b28 * PDM[3][1] / b32).ln();
                let hc32 = PDM[3][5] * PDL[1][7];
                let zc32 = PDM[3][4] * PDL[1][6];
                output.molecular_oxygen *= ccor(z, rl, hc32, zc32);
            }
            // Correction for a general departure from diffusive equilibrium above the
            // lower boundary, which applies at every altitude.
            let hcc32 = PDM[3][7] * PDL[1][22];
            let hcc232 = PDM[3][7] * PDL[0][22];
            let zcc32 = PDM[3][6] * PDL[1][21];
            let rc32 = ((input.f107_average - 150.0) * (sw[1] * PDL[0][23]) + 1.0)
                * (PDM[3][3] * PDL[1][23]);
            output.molecular_oxygen *= ccor2(z, rc32, hcc32, zcc32, hcc232);
        }

        // Argon.
        let g40 = globe7(&PD[5], &self.phases.pd[5], input, &self.resolved, work) * sw[21];
        let db40 = g40.exp() * PDM[4][0] * PD[5][0];
        output.argon = densu(
            work, z, db40, tinf, tlb, 40.0, ALPHA[4], &mut tz, PTM[5], s, &zn1,
        );
        if sw[15] != 0.0 && z.re() <= ALTL[4] {
            let zh40 = PDM[4][2];
            let b40 = densu(
                work,
                T::lift(zh40),
                db40,
                tinf,
                tlb,
                40.0 - xmm,
                ALPHA[4] - 1.0,
                &mut tz,
                PTM[5],
                s,
                &zn1,
            );
            let dm40 = densu(work, z, b40, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
            work.dm40 = dm40;
            let zhm40 = zhm28;
            output.argon = dnet(output.argon, dm40, zhm40, xmm, 40.0);
            let rl = (b28 * PDM[4][1] / b40).ln();
            let hc40 = PDM[4][5] * PDL[1][9];
            let zc40 = PDM[4][4] * PDL[1][8];
            output.argon *= ccor(z, rl, hc40, zc40);
        }

        // Hydrogen.
        let g1 = globe7(&PD[6], &self.phases.pd[6], input, &self.resolved, work) * sw[21];
        let db01 = g1.exp() * PDM[5][0] * PD[6][0];
        output.hydrogen = densu(
            work, z, db01, tinf, tlb, 1.0, ALPHA[6], &mut tz, PTM[5], s, &zn1,
        );
        if sw[15] != 0.0 && z.re() <= ALTL[6] {
            let zh01 = PDM[5][2];
            let b01 = densu(
                work,
                T::lift(zh01),
                db01,
                tinf,
                tlb,
                1.0 - xmm,
                ALPHA[6] - 1.0,
                &mut tz,
                PTM[5],
                s,
                &zn1,
            );
            let dm01 = densu(work, z, b01, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
            work.dm01 = dm01;
            let zhm01 = zhm28;
            output.hydrogen = dnet(output.hydrogen, dm01, zhm01, xmm, 1.0);
            // The reference writes sqrt(pdl[1][17]*pdl[1][17]) for the magnitude.
            let rl = (b28 * PDM[5][1] * PDL[1][17].abs() / b01).ln();
            let hc01 = PDM[5][5] * PDL[1][11];
            let zc01 = PDM[5][4] * PDL[1][10];
            output.hydrogen *= ccor(z, rl, hc01, zc01);
            // Chemistry correction.
            let hcc01 = PDM[5][7] * PDL[1][19];
            let zcc01 = PDM[5][6] * PDL[1][18];
            let rc01 = PDM[5][3] * PDL[1][20];
            output.hydrogen *= ccor(z, T::lift(rc01), hcc01, zcc01);
        }

        // Atomic nitrogen.
        let g14 = globe7(&PD[7], &self.phases.pd[7], input, &self.resolved, work) * sw[21];
        let db14 = g14.exp() * PDM[6][0] * PD[7][0];
        output.atomic_nitrogen = densu(
            work, z, db14, tinf, tlb, 14.0, ALPHA[7], &mut tz, PTM[5], s, &zn1,
        );
        if sw[15] != 0.0 && z.re() <= ALTL[7] {
            let zh14 = PDM[6][2];
            let b14 = densu(
                work,
                T::lift(zh14),
                db14,
                tinf,
                tlb,
                14.0 - xmm,
                ALPHA[7] - 1.0,
                &mut tz,
                PTM[5],
                s,
                &zn1,
            );
            let dm14 = densu(work, z, b14, tinf, tlb, xmm, 0.0, &mut tz, PTM[5], s, &zn1);
            work.dm14 = dm14;
            let zhm14 = zhm28;
            output.atomic_nitrogen = dnet(output.atomic_nitrogen, dm14, zhm14, xmm, 14.0);
            let rl = (b28 * PDM[6][1] * PDL[0][2].abs() / b14).ln();
            let hc14 = PDM[6][5] * PDL[0][1];
            let zc14 = PDM[6][4] * PDL[0][0];
            output.atomic_nitrogen *= ccor(z, rl, hc14, zc14);
            // Chemistry correction.
            let hcc14 = PDM[6][7] * PDL[0][4];
            let zcc14 = PDM[6][6] * PDL[0][3];
            let rc14 = PDM[6][3] * PDL[0][5];
            output.atomic_nitrogen *= ccor(z, T::lift(rc14), hcc14, zcc14);
        }

        // Anomalous oxygen: hot atomic oxygen, on its own barometric profile rather than
        // the diffusive one, and with no mixed correction.
        let g16h = globe7(&PD[8], &self.phases.pd[8], input, &self.resolved, work) * sw[21];
        let db16h = g16h.exp() * PDM[7][0] * PD[8][0];
        let tho = PDM[7][9] * PDL[0][6];
        let dd = densu(
            work,
            z,
            db16h,
            T::lift(tho),
            T::lift(tho),
            16.0,
            ALPHA[8],
            &mut tz,
            PTM[5],
            s,
            &zn1,
        );
        let zsht = PDM[7][5];
        let zmho = PDM[7][4];
        let zsho = scalh(work, zmho, 16.0, T::lift(tho));
        output.anomalous_oxygen =
            dd * ((T::lift(-zsht) / zsho) * ((-(z - zmho) / zsht).exp() - 1.0)).exp();

        output.total_mass = mass_density(&output, false);

        // The temperature at the altitude, from the same profile with no mass.
        let z = alt.abs();
        densu(
            work,
            z,
            T::lift(1.0),
            tinf,
            tlb,
            0.0,
            0.0,
            &mut tz,
            PTM[5],
            s,
            &zn1,
        );
        output.temperature = tz;

        if self.switches.units == Units::Si {
            to_si(&mut output);
        }
        output
    }
}

/// The mass density summed over the species, in the model's own units.
///
/// Which species it sums is the only difference between `gtd7` and `gtd7d`: the latter
/// includes anomalous oxygen.
fn mass_density<T: Scalar>(output: &Output<T>, anomalous: bool) -> T {
    let sum = output.helium * 4.0
        + output.atomic_oxygen * 16.0
        + output.molecular_nitrogen * 28.0
        + output.molecular_oxygen * 32.0
        + output.argon * 40.0
        + output.hydrogen
        + output.atomic_nitrogen * 14.0;
    let sum = if anomalous {
        sum + output.anomalous_oxygen * 16.0
    } else {
        sum
    };
    sum * NUCLEON_MASS
}

/// The reference's switch 0, applied where it applies it.
///
/// Every density, the mass density among them, is scaled from cm^-3 to m^-3, and then the
/// mass density is divided by a thousand. Both operations, in that order: g/cm^3 scaled by
/// 1e6 and then by 1e-3 is kg/m^3, and doing it in one step would not round the same way.
fn to_si<T: Scalar>(output: &mut Output<T>) {
    output.helium *= 1.0e6;
    output.atomic_oxygen *= 1.0e6;
    output.molecular_nitrogen *= 1.0e6;
    output.molecular_oxygen *= 1.0e6;
    output.argon *= 1.0e6;
    output.total_mass *= 1.0e6;
    output.hydrogen *= 1.0e6;
    output.atomic_nitrogen *= 1.0e6;
    output.anomalous_oxygen *= 1.0e6;
    output.total_mass /= 1000.0;
}

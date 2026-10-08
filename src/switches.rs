//! Which variations the model includes, and which units it answers in.
//!
//! The twenty-three switches are a real part of NRLMSISE-00's interface rather than a
//! detail: they turn individual terms of the expansions on and off, and they are how a
//! published result that used a reduced model can be reproduced. They are named here
//! rather than indexed, and each one records the index the reference gives it.
//!
//! A switch has three states, not two. [`Switch::CrossTermsOnly`] is the reference's
//! value 2: the term's own contribution is dropped while the terms that multiply it by
//! something else are kept.
//!
//! # Magnetic activity is not here
//!
//! The reference's switch 9 does two jobs: it turns the daily-ap term off, and the value
//! -1 selects the seven-element ap history instead of a single daily value. The second is
//! a statement about the data, not about the model, so it lives on the input as
//! [`crate::MagneticActivity`] and the data cannot disagree with the mode.

/// Whether one variation is included.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Switch {
    /// Left out, along with every cross term that involves it. The reference's 0.
    Off,
    /// Included. The reference's 1.
    #[default]
    On,
    /// The term itself left out, but cross terms involving it kept. The reference's 2.
    CrossTermsOnly,
}

impl Switch {
    /// The reference's `sw`: one when the term itself is included.
    fn sw(self) -> f64 {
        f64::from(self == Switch::On)
    }

    /// The reference's `swc`: one when cross terms involving it are included.
    fn swc(self) -> f64 {
        f64::from(self != Switch::Off)
    }
}

/// What units the densities and the mass density come out in.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Units {
    /// Number densities in m^-3 and mass density in kg/m^3.
    ///
    /// The reference's switch 0 set, and this crate's default. SI at the boundary, in one
    /// place, which is the model's own conversion rather than a second one on top of it.
    #[default]
    Si,
    /// Number densities in cm^-3 and mass density in g/cm^3.
    ///
    /// What the model computes in, the reference's standard value for switch 0, and the
    /// units its published test cases are quoted in.
    Native,
}

/// The twenty-three variations, and the units.
///
/// [`Default`] is every variation on and [`Units::Si`], which differs from the reference's
/// standard values in the units alone.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Switches {
    /// Switch 0 in the reference, where it selects meters and kilograms.
    pub units: Units,
    /// 1: the F10.7 effect on the mean.
    pub f107_mean: Switch,
    /// 2: the time independent terms. Also drops the latitude variation of gravity.
    pub time_independent: Switch,
    /// 3: the symmetrical annual variation.
    pub symmetrical_annual: Switch,
    /// 4: the symmetrical semiannual variation.
    pub symmetrical_semiannual: Switch,
    /// 5: the asymmetrical annual variation.
    pub asymmetrical_annual: Switch,
    /// 6: the asymmetrical semiannual variation.
    pub asymmetrical_semiannual: Switch,
    /// 7: the diurnal variation.
    pub diurnal: Switch,
    /// 8: the semidiurnal variation.
    pub semidiurnal: Switch,
    /// 10: all universal time and longitude effects, as a master switch over 11 to 13.
    pub ut_and_longitude: Switch,
    /// 11: the longitudinal variation.
    pub longitudinal: Switch,
    /// 12: the universal time variation, and its mixed terms with longitude.
    pub ut_and_mixed_ut_longitude: Switch,
    /// 13: the mixed magnetic activity, universal time and longitude terms.
    pub mixed_ap_ut_longitude: Switch,
    /// 14: the terdiurnal variation.
    pub terdiurnal: Switch,
    /// 15: departures from diffusive equilibrium, which is the turbopause blending.
    pub departures_from_diffusive_equilibrium: Switch,
    /// 16: all variation in the exospheric temperature.
    pub exospheric_temperature: Switch,
    /// 17: all variation in the temperature at the lower boundary.
    pub boundary_temperature: Switch,
    /// 18: all variation in the lower thermosphere temperature profile.
    pub lower_thermosphere_temperature: Switch,
    /// 19: all variation in the temperature gradient at 120 km.
    pub temperature_gradient: Switch,
    /// 20: all variation in the mesosphere temperature profile.
    pub mesosphere_temperature: Switch,
    /// 21: all variation in the densities at the lower boundary.
    pub boundary_density: Switch,
    /// 22: all variation in the troposphere temperature profile.
    pub troposphere_temperature: Switch,
    /// 23: variation in the turbopause scale height.
    pub turbopause_scale_height: Switch,
}

impl Switches {
    /// The reference's standard values: every variation on, in the model's own units.
    ///
    /// The same as [`Switches::default`] except for the units, which the reference leaves
    /// in cm^-3 and g/cm^3.
    pub fn reference() -> Self {
        Self {
            units: Units::Native,
            ..Self::default()
        }
    }

    /// Everything off, as a base to turn individual variations back on from.
    pub fn none() -> Self {
        Self {
            units: Units::default(),
            f107_mean: Switch::Off,
            time_independent: Switch::Off,
            symmetrical_annual: Switch::Off,
            symmetrical_semiannual: Switch::Off,
            asymmetrical_annual: Switch::Off,
            asymmetrical_semiannual: Switch::Off,
            diurnal: Switch::Off,
            semidiurnal: Switch::Off,
            ut_and_longitude: Switch::Off,
            longitudinal: Switch::Off,
            ut_and_mixed_ut_longitude: Switch::Off,
            mixed_ap_ut_longitude: Switch::Off,
            terdiurnal: Switch::Off,
            departures_from_diffusive_equilibrium: Switch::Off,
            exospheric_temperature: Switch::Off,
            boundary_temperature: Switch::Off,
            lower_thermosphere_temperature: Switch::Off,
            temperature_gradient: Switch::Off,
            mesosphere_temperature: Switch::Off,
            boundary_density: Switch::Off,
            troposphere_temperature: Switch::Off,
            turbopause_scale_height: Switch::Off,
        }
    }

    /// The switches in the reference's own order, so they can be resolved by index.
    ///
    /// Index 9 is absent from this list; see the module docs.
    fn in_order(&self) -> [(usize, Switch); 22] {
        [
            (1, self.f107_mean),
            (2, self.time_independent),
            (3, self.symmetrical_annual),
            (4, self.symmetrical_semiannual),
            (5, self.asymmetrical_annual),
            (6, self.asymmetrical_semiannual),
            (7, self.diurnal),
            (8, self.semidiurnal),
            (10, self.ut_and_longitude),
            (11, self.longitudinal),
            (12, self.ut_and_mixed_ut_longitude),
            (13, self.mixed_ap_ut_longitude),
            (14, self.terdiurnal),
            (15, self.departures_from_diffusive_equilibrium),
            (16, self.exospheric_temperature),
            (17, self.boundary_temperature),
            (18, self.lower_thermosphere_temperature),
            (19, self.temperature_gradient),
            (20, self.mesosphere_temperature),
            (21, self.boundary_density),
            (22, self.troposphere_temperature),
            (23, self.turbopause_scale_height),
        ]
    }

    /// The switches as the reference's `sw` and `swc` arrays.
    ///
    /// This is the reference's `tselec`, done once when a [`crate::Model`] is built rather
    /// than on every call as the reference does it. Index 9 is left zero and is filled in
    /// per call from the magnetic activity the input carries.
    pub(crate) fn resolve(&self) -> Resolved {
        let mut sw = [0.0; 24];
        let mut swc = [0.0; 24];
        sw[0] = f64::from(self.units == Units::Si);
        swc[0] = sw[0];
        for (index, switch) in self.in_order() {
            sw[index] = switch.sw();
            swc[index] = switch.swc();
        }
        Resolved { sw, swc }
    }
}

/// The switches in the form the expansions multiply by.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Resolved {
    /// Whether each term itself is included, as 1.0 or 0.0.
    pub sw: [f64; 24],
    /// Whether cross terms involving it are included, as 1.0 or 0.0.
    pub swc: [f64; 24],
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three states resolve the way the reference's `tselec` resolves them.
    #[test]
    fn test_the_switch_states_resolve_as_tselec_does() {
        assert_eq!((Switch::Off.sw(), Switch::Off.swc()), (0.0, 0.0));
        assert_eq!((Switch::On.sw(), Switch::On.swc()), (1.0, 1.0));
        // The reference's 2: sw is zero because it is not 1, swc is one because it is > 0.
        assert_eq!(
            (Switch::CrossTermsOnly.sw(), Switch::CrossTermsOnly.swc()),
            (0.0, 1.0)
        );
    }

    /// Every index but 9 is covered exactly once.
    #[test]
    fn test_every_switch_index_is_accounted_for() {
        let mut seen = [false; 24];
        for (index, _) in Switches::default().in_order() {
            assert!(!seen[index], "index {index} twice");
            seen[index] = true;
        }
        for (index, covered) in seen.iter().enumerate() {
            let expected = index != 0 && index != 9;
            assert_eq!(*covered, expected, "index {index}");
        }
    }

    #[test]
    fn test_the_default_is_everything_on_in_si() {
        let resolved = Switches::default().resolve();
        assert_eq!(resolved.sw[0], 1.0);
        for index in 1..24 {
            if index == 9 {
                continue;
            }
            assert_eq!(resolved.sw[index], 1.0, "sw[{index}]");
            assert_eq!(resolved.swc[index], 1.0, "swc[{index}]");
        }
        // The reference's own standard values differ only in the units.
        assert_eq!(Switches::reference().resolve().sw[0], 0.0);
    }
}

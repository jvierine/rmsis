//! What the model is asked about, and what it answers with.

use crate::scalar::Scalar;

/// The geomagnetic activity the model reads, and which of its two modes that selects.
///
/// The reference carries a daily value and a seven-element history side by side and picks
/// between them with switch 9, which means the switch and the data can disagree. Here the
/// data is the mode.
///
/// [`MagneticActivity::History`] is what a real propagation uses: a storm's effect on
/// density lags it by hours, which a single daily mean cannot express. The reference's
/// published cases 15 and 16 are the two that exercise it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MagneticActivity<T> {
    /// The daily-ap term left out, which is the reference's switch 9 set to zero.
    ///
    /// The daily value is still carried because the reference still reads it for the mixed
    /// magnetic, universal time and longitude term, which switch 13 governs separately.
    Off {
        /// The daily mean ap.
        daily_ap: T,
    },
    /// The daily mean ap alone, which is the reference's switch 9 set to one.
    Daily {
        /// The daily mean ap.
        ap: T,
    },
    /// The seven-element history, which is the reference's switch 9 set to -1.
    History {
        /// In the order the model wants them:
        ///
        /// 0. the daily mean ap
        /// 1. the three-hour ap for the interval the epoch falls in
        /// 2. three hours before that
        /// 3. six hours before
        /// 4. nine hours before
        /// 5. the mean of the eight intervals from 12 to 33 hours before
        /// 6. the mean of the eight intervals from 36 to 57 hours before
        ap: [T; 7],
    },
}

impl<T: Scalar> MagneticActivity<T> {
    /// The daily mean, which every variant carries.
    pub(crate) fn daily(&self) -> T {
        match self {
            Self::Off { daily_ap } => *daily_ap,
            Self::Daily { ap } => *ap,
            Self::History { ap } => ap[0],
        }
    }

    /// The reference's `sw[9]`: 0 off, 1 daily, -1 history.
    pub(crate) fn sw(&self) -> f64 {
        match self {
            Self::Off { .. } => 0.0,
            Self::Daily { .. } => 1.0,
            Self::History { .. } => -1.0,
        }
    }
}

/// Where and when the model is being asked about.
///
/// Altitude is geodetic and in km, latitude and longitude are geodetic and in degrees, and
/// the time of day is in seconds of universal time. These are the reference's own units at
/// the input; only the output is converted.
///
/// # Local solar time
///
/// Universal time, longitude and local solar time enter the model independently, and it
/// does not check that they agree. They should: `local_solar_time_h = seconds_of_day/3600 +
/// longitude_deg/15`. [`Input::consistent`] works it out so it cannot be forgotten. The
/// field is left public because the reference's published cases deliberately break the
/// relation, and reproducing them is the point of the fixtures.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Input<T> {
    /// Day of the year, 1 to 366.
    ///
    /// An integer, as in the reference, so it carries no derivative. Nothing in the model
    /// differentiates with respect to the day.
    pub day_of_year: i32,
    /// Seconds of the day, universal time.
    pub seconds_of_day: T,
    /// Geodetic altitude, km.
    pub altitude_km: T,
    /// Geodetic latitude, degrees.
    pub latitude_deg: T,
    /// Geodetic longitude, degrees east.
    pub longitude_deg: T,
    /// Local apparent solar time, hours. See the note on [`Input`].
    pub local_solar_time_h: T,
    /// The 81-day average of the F10.7 solar flux, centered on the day.
    ///
    /// At the flux's value at the Earth's actual distance from the Sun, not at 1 AU.
    pub f107_average: T,
    /// The daily F10.7 solar flux of the *previous* day.
    pub f107: T,
    /// The geomagnetic activity, and which mode it selects.
    pub magnetic: MagneticActivity<T>,
}

impl<T: Scalar> Input<T> {
    /// An input whose local solar time agrees with its universal time and longitude.
    ///
    /// `local_solar_time_h = seconds_of_day/3600 + longitude_deg/15`, which is what the
    /// model's own documentation asks for.
    #[allow(clippy::too_many_arguments)]
    pub fn consistent(
        day_of_year: i32,
        seconds_of_day: T,
        altitude_km: T,
        latitude_deg: T,
        longitude_deg: T,
        f107_average: T,
        f107: T,
        magnetic: MagneticActivity<T>,
    ) -> Self {
        Self {
            day_of_year,
            seconds_of_day,
            altitude_km,
            latitude_deg,
            longitude_deg,
            local_solar_time_h: seconds_of_day / 3600.0 + longitude_deg / 15.0,
            f107_average,
            f107,
            magnetic,
        }
    }
}

/// What the model answers: eight species, a total, and two temperatures.
///
/// The number densities are in m^-3 and the mass density in kg/m^3 under
/// [`crate::Units::Si`], or cm^-3 and g/cm^3 under [`crate::Units::Native`].
///
/// Atomic oxygen, hydrogen and nitrogen are zero below 72.5 km, where the model does not
/// carry them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Output<T> {
    /// Helium number density.
    pub helium: T,
    /// Atomic oxygen number density.
    pub atomic_oxygen: T,
    /// Molecular nitrogen number density.
    pub molecular_nitrogen: T,
    /// Molecular oxygen number density.
    pub molecular_oxygen: T,
    /// Argon number density.
    pub argon: T,
    /// Total mass density.
    ///
    /// Which species this sums over is the difference between [`crate::Model::gtd7`] and
    /// [`crate::Model::gtd7d`], and it is the one number a drag calculation wants. See
    /// those two.
    pub total_mass: T,
    /// Atomic hydrogen number density.
    pub hydrogen: T,
    /// Atomic nitrogen number density.
    pub atomic_nitrogen: T,
    /// Anomalous oxygen number density.
    ///
    /// Hot atomic oxygen, which matters for drag above about 500 km and is the species
    /// `gtd7` leaves out of the total while `gtd7d` includes it.
    pub anomalous_oxygen: T,
    /// Exospheric temperature, K.
    ///
    /// Set to the global average below 120 km, where the model does not vary it.
    pub exospheric_temperature: T,
    /// Temperature at the altitude asked for, K.
    pub temperature: T,
}

impl<T: Scalar> Output<T> {
    /// All zero.
    pub(crate) fn zeroed() -> Self {
        let zero = T::lift(0.0);
        Self {
            helium: zero,
            atomic_oxygen: zero,
            molecular_nitrogen: zero,
            molecular_oxygen: zero,
            argon: zero,
            total_mass: zero,
            hydrogen: zero,
            atomic_nitrogen: zero,
            anomalous_oxygen: zero,
            exospheric_temperature: zero,
            temperature: zero,
        }
    }

    /// The nine densities in the order the reference's `d` array holds them.
    ///
    /// He, O, N2, O2, Ar, total mass, H, N, anomalous O. The order the published test
    /// cases are quoted in, which is why it is worth having.
    pub fn densities(&self) -> [T; 9] {
        [
            self.helium,
            self.atomic_oxygen,
            self.molecular_nitrogen,
            self.molecular_oxygen,
            self.argon,
            self.total_mass,
            self.hydrogen,
            self.atomic_nitrogen,
            self.anomalous_oxygen,
        ]
    }
}

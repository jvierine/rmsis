//! Optional Rust interface to an unmodified, separately installed NRLMSIS 2.1 library.
//! Build the C ABI adapter with `scripts/build_msis21.py`; NRL sources and
//! parameters are not bundled. The original NRL license applies to that backend.
use libloading::Library;
use std::{ffi::CString, path::Path, sync::Mutex};

type Init = unsafe extern "C" fn(*const std::ffi::c_char, i32, i32);
type Calc = unsafe extern "C" fn(*const f64, *const f64, *mut f64);
static STATE: Mutex<Option<(usize, String, bool)>> = Mutex::new(None);

/// Safe scalar input; Ap order follows pymsis/NRL's seven-element convention.
pub struct Input {
    /// UTC day of year, 1 through 366.
    pub day_of_year: f64,
    /// UTC seconds after midnight.
    pub seconds_of_day: f64,
    /// Geodetic altitude, km.
    pub altitude_km: f64,
    /// Geodetic latitude, degrees north.
    pub latitude_deg: f64,
    /// Longitude, degrees east.
    pub longitude_deg: f64,
    /// 81-day mean F10.7 solar flux, sfu.
    pub f107_average: f64,
    /// Previous day's F10.7 solar flux, sfu.
    pub f107: f64,
    /// Daily, current/previous three 3-hour Ap, then 12–33h and 36–57h means.
    pub ap: [f64; 7],
    /// Use the seven-element storm-time history instead of daily Ap alone.
    pub storm_time: bool,
}
/// SI densities in pymsis order and temperatures in kelvin.
pub struct Output {
    /// rho kg/m^3; N2,O2,O,He,H,Ar,N,anomalous O,NO in m^-3.
    /// Undefined species are None (the original library's missing-value sentinel).
    pub densities: [Option<f64>; 10],
    /// Local neutral temperature, K.
    pub temperature_k: f64,
    /// Exospheric temperature, K.
    pub exospheric_temperature_k: f64,
}
/// A loaded backend. Calls are serialized because the original model has saved state.
pub struct Model {
    _library: Library,
    init: Init,
    calc: Calc,
    parameters: CString,
}
impl Model {
    /// Load the compiled adapter and an existing original `msis21.parm` file.
    ///
    /// # Safety
    /// The shared library must be trusted and export the adapter ABI provided in
    /// `ffi/msis21_bridge.f90`. Do not load arbitrary third-party binaries.
    pub unsafe fn load(
        library: impl AsRef<Path>,
        parameters: impl AsRef<Path>,
    ) -> Result<Self, String> {
        let parameter_path = parameters
            .as_ref()
            .canonicalize()
            .map_err(|e| e.to_string())?;
        let text = parameter_path
            .to_str()
            .ok_or("Parameter path must be UTF-8")?;
        if text.len() > 128 {
            return Err("NRL limits the parameter filename to 128 bytes".into());
        }
        if !parameter_path.is_file() {
            return Err("Parameter file is not a regular file".into());
        }
        let parameters = CString::new(text).map_err(|e| e.to_string())?;
        let loaded = Library::new(library.as_ref()).map_err(|e| e.to_string())?;
        let init = *loaded
            .get::<Init>(b"rmsis21_init\0")
            .map_err(|e| e.to_string())?;
        let calc = *loaded
            .get::<Calc>(b"rmsis21_calc\0")
            .map_err(|e| e.to_string())?;
        Ok(Self {
            _library: loaded,
            init,
            calc,
            parameters,
        })
    }
    /// Evaluate the original model without calling Python or accessing the network.
    pub fn evaluate(&self, input: &Input) -> Result<Output, String> {
        let values = [
            input.day_of_year,
            input.seconds_of_day,
            input.altitude_km,
            input.latitude_deg,
            input.longitude_deg,
            input.f107_average,
            input.f107,
        ];
        if !values.iter().chain(input.ap.iter()).all(|x| x.is_finite())
            || !(1. ..367.).contains(&input.day_of_year)
            || !(0. ..86400.).contains(&input.seconds_of_day)
            || !(0. ..=2000.).contains(&input.altitude_km)
            || !(-90. ..=90.).contains(&input.latitude_deg)
            || input.f107 <= 0.
            || input.f107_average <= 0.
            || input.ap.iter().any(|x| *x < 0.)
        {
            return Err("Invalid NRLMSIS input".into());
        }
        let mut state = STATE.lock().map_err(|_| "Backend lock poisoned")?;
        let key = (
            self.calc as usize,
            self.parameters.to_string_lossy().into_owned(),
            input.storm_time,
        );
        let mut output = [0.; 12];
        // SAFETY: validated inputs, exact adapter signatures, held library and
        // process-wide mutex protecting both initialization and evaluation.
        unsafe {
            if state.as_ref() != Some(&key) {
                (self.init)(
                    self.parameters.as_ptr(),
                    self.parameters.as_bytes().len() as i32,
                    i32::from(input.storm_time),
                );
                *state = Some(key);
            }
            (self.calc)(values.as_ptr(), input.ap.as_ptr(), output.as_mut_ptr());
        }
        Ok(Output {
            densities: std::array::from_fn(|i| {
                if output[i] < 0. {
                    None
                } else {
                    Some(output[i])
                }
            }),
            temperature_k: output[10],
            exospheric_temperature_k: output[11],
        })
    }
}
impl Drop for Model {
    fn drop(&mut self) {
        // Clear the cached identity before unloading, since a future loader may
        // reuse the same address for a different instance.
        if let Ok(mut state) = STATE.lock() {
            *state = None;
        }
    }
}

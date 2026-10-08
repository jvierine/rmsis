# rmsis

Rust atmospheric-model APIs: **NRLMSISE-00** in pure Rust (default), and an
optional **NRLMSIS 2.1 interface linked to the original Fortran library**.
Neither backend calls Python, Git, or the network. Drivers (F10.7 and Ap) are
explicit inputs. Only the optional 2.1 backend requires Fortran.

```rust
use rmsis::{Input, MagneticActivity, Model};
let input = Input::consistent(172, 43200.0, 100.0, 69.0, 19.0,
    80.0, 80.0, MagneticActivity::Daily { ap: 4.0 });
let output = Model::new().gtd7d(&input);
println!("rho={} kg/m^3; T={} K", output.total_mass, output.temperature);
```

Coordinates: geodetic latitude/longitude in degrees, altitude in km, universal
time in seconds, day of year 1–366. `Input::consistent` derives local solar
time from UT and east longitude. Default outputs are SI: number densities
m^-3, mass density kg/m^3, temperature K. `MagneticActivity::History` accepts
the standard seven Ap values and selects the storm-time mode. `gtd7d` includes
anomalous oxygen in drag density; `gtd7` excludes it. Individual species and
both local and exospheric temperatures are available in `Output`.

## Model version and origin

The default engine implements **NRLMSISE-00**, compared to `pymsis` with
`version=0`; it has no nitric oxide output. The `msis21` feature exposes the
original 2.1 model, including NO, through a separate C-ABI adapter; it is not
a translation of 2.1. No NRL 2.1 source, parameters, or binaries are bundled.
Rust code is MIT; the separately obtained NRL backend retains its own license.

The complete engine is adapted from Connor Johnstone's MIT-licensed
[`nrlmsise00` 0.2.0](https://crates.io/crates/nrlmsise00/0.2.0), retaining its
generic scalar support, switches, and coefficient tables. This is an attributed
adaptation, not an independently authored atmospheric model. Upstream license
is in `licenses/nrlmsise00-MIT.txt`. The helium correction at exactly 200 km
includes the boundary, matching the original Fortran rather than the upstream
C-derived strict inequality. Physical model: Picone, Hedin, Drob, and
Aikin (2002), DOI 10.1029/2002JA009430; original public-domain C implementation
by Dominik Brodowski. Optimizations under `fast` are optional and must satisfy
the same reference tests.

## Verification

```sh
cargo test --all-features
conda run -n base python -m pip install pymsis==0.11.0 numpy h5py
conda run -n base python scripts/generate_reference.py
cargo test --all-features
```

Python is only needed to regenerate reference fixtures. Checked-in JSON test
fixtures let ordinary `cargo test` run without Python; HDF5 stores the numerical
reference product and provenance. 430 cases per backend span seasons, UT, latitudes, altitudes,
solar activity, and daily/history Ap. Tolerances account for pymsis's float32
wrapper: MSISE-00 uses 0.015% relative density tolerance and 0.01 K temperature
tolerance, plus negligible absolute floors. The comparison does not establish
empirical accuracy of the model. Ordinary tests skip the externally linked 2.1
test; it was separately run against an unmodified backend before publication.

## Link the original NRLMSIS 2.1 library

Obtain the [NRL distribution](https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/)
and review its license before using it. Install `gfortran`, then:

```sh
conda run -n base python scripts/build_msis21.py /path/to/NRLMSIS2.1
```

The script compiles the original sources without edits, plus our separate
`ffi/msis21_bridge.f90`. It prints the library path (`.so`, `.dylib`, or `.dll`).
Enable `msis21` in Cargo and load that path and the original `msis21.parm`:

```rust,no_run
# #[cfg(feature = "msis21")]
# fn main() -> Result<(), String> {
use rmsis::msis21::{Input, Model};
// SAFETY: load only a trusted library built with the supplied adapter.
let model = unsafe { Model::load("/path/to/librmsis21.so", "/path/to/msis21.parm") }?;
let input = Input { day_of_year: 172., seconds_of_day: 43200., altitude_km: 100.,
    latitude_deg: 69., longitude_deg: 19., f107_average: 80., f107: 80.,
    ap: [4.; 7], storm_time: false };
let output = model.evaluate(&input)?;
println!("T={} K; rho={:?} kg/m^3", output.temperature_k, output.densities[0]);
# Ok(()) }
# #[cfg(not(feature = "msis21"))] fn main() {}
```

All model calls are serialized to protect the original saved state, including
switch changes between daily and storm-time Ap. This lock cannot protect calls
made directly to the same Fortran library by unrelated code. Parameter paths
are limited to 128 UTF-8 bytes by the original model. The adapter recomputes the
mass-density weighted sum from returned species to avoid the original routine's
uninitialized self-term; original NRL source remains untouched.

```sh
RMSIS21_LIBRARY=/absolute/path/to/librmsis21.so \
RMSIS21_PARAMETERS=/absolute/path/to/msis21.parm \
cargo test --features msis21 --test msis21_reference -- --ignored
```

The 2.1 comparison uses 0.005% relative density and 0.003 K temperature tolerance.
Library loading is explicitly unsafe because it executes a supplied binary;
validated model evaluation is safe Rust. On Windows use a matching-architecture
GNU Fortran toolchain and make its runtime DLLs available on `PATH`.

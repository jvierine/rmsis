# Validation record — 2026-10-08

On macOS arm64, using Cargo 1.97.1, GNU Fortran, and pymsis 0.11.0:

- `cargo test --all-features`: unit tests, 430 NRLMSISE-00 reference cases,
  and documentation examples passed.
- The explicitly enabled linked NRLMSIS 2.1 test passed all 430 reference
  cases, alternating daily/history Ap and two Rust model handles.
- `cargo clippy --all-targets --all-features` passed without warnings.

The linked test used the official NRL archive, downloaded from
https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/nrlmsis2.1.tar.gz
and compiled without editing its files. Original source and parameter data
are not included in this repository. SHA-256:

```
41e47b29f795d36a5cc252b2858aa2a384c4a7323ace3d48d3ea2f2b37a1a6a8  nrlmsis2.1.tar.gz
a322a749f368e73117dd20f3fdcf7389dabc5509f4c27073cc5580999381b508  msis21.parm
```

The parameter-file hash also matched pymsis's installed parameter file.
`scripts/generate_reference.py` explicitly supplies space-weather inputs;
reference generation does not download drivers. HDF5 products preserve the
reference values, inputs, seed, model identifier, and Python-package version.
JSON mirrors are included only so Rust tests need no Python/HDF5 dependency.

Tolerances and the adapter-level mass-density safeguard are documented in
README.md. Numerical agreement is not validation of physical model accuracy.
The CI matrix checks ordinary Rust tests on Linux/macOS/Windows; external
Fortran tests require separately supplied original distributions and are
not claimed as having been run on Windows or Linux in this validation record.

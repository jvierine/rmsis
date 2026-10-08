#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

#[cfg(feature = "msis21")]
pub mod msis21;
pub mod scalar;

mod expansion;
mod helpers;
mod input;
mod model;
mod profile;
mod switches;
mod tables;
mod work;

pub use input::{Input, MagneticActivity, Output};
pub use model::Model;
pub use scalar::Scalar;
pub use switches::{Switch, Switches, Units};

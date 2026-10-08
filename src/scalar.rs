//! The scalar type the model is evaluated in.
//!
//! Instantiating the model at [`f64`] gives density. Instantiating it at a forward-mode
//! dual number gives density *and* its exact derivatives, which is the reason this crate
//! exists: density falls by a factor of e about every 50 km, so a drag partial that drops
//! its dependence on position is badly wrong rather than slightly wrong.
//!
//! # It is `num-dual`'s trait, not one of ours
//!
//! [`Scalar`] declares one method of its own and inherits the rest from
//! [`num_dual::DualNum`], which `num-dual` implements for `f64` and `f32` as well as for
//! its dual types. It is then blanket-implemented, so **every** type satisfying that trait
//! at `f64` precision is a scalar here, including nestings and types this crate has never
//! heard of.
//!
//! That is deliberate. An earlier version declared `exp`, `ln`, `powf` and the rest itself
//! and implemented them for a list of types. Two things were wrong with that. The list
//! could not be complete, and nobody else could extend it: a caller whose scalar was a
//! dual type the list had missed was stopped by the orphan rule, with no way to implement
//! this crate's trait for it. And declaring the same methods a second time made this trait
//! incompatible with every other crate's trait over the same types — `ordinary-diffeq`
//! declares one of the same shape — so a caller generic over both was adapting between
//! two names for one thing. Taking the bound from `num-dual` leaves one name.
//!
//! A caller generic over another crate's scalar trait still writes both bounds, since
//! neither trait implies the other, but the second one is now a formality that any dual
//! number satisfies rather than a wall that some do not.
//!
//! # What has to be in scope
//!
//! Nothing, to write generic code. The methods belong to [`Scalar`]'s supertraits, and a
//! type parameter bounded by `Scalar` carries them without the supertraits being imported:
//!
//! ```
//! use rmsis::Scalar;
//!
//! fn scale_height<T: Scalar>(density: T, slope: T) -> T {
//!     // ln from DualNum, abs from Signed, re from DualStruct, lift from Scalar, and only
//!     // Scalar is imported.
//!     let _ = (density.ln(), density.abs(), density.re());
//!     -density / slope + T::lift(0.0)
//! }
//!
//! assert!((scale_height(1.0_f64, -0.02) - 50.0).abs() < 1e-12);
//! ```
//!
//! Calling the same methods on a *concrete* dual number does need the trait that declares
//! them, and the three are re-exported here for that: [`DualNum`] for the arithmetic,
//! [`DualStruct`] for [`DualStruct::re`], [`Signed`] for [`Signed::abs`].
//!
//! # A literal-only call has to name its scalar
//!
//! The one wrinkle of the bound being a blanket one. Inference will not settle a float
//! literal on `f64` when the only thing constraining it is `T: Scalar`, because the
//! blanket implementation admits every dual number too, so a call built entirely from
//! literals is ambiguous and says so at the first method call on its result. Name the
//! scalar once and the whole call follows:
//!
//! ```
//! # use rmsis::{Input, MagneticActivity, Model};
//! let model = Model::new();
//! let input = Input::consistent(
//!     172, 29000.0, 400.0, 60.0, -70.0, 150.0, 150.0,
//!     MagneticActivity::Daily { ap: 4.0 },
//! );
//! let output = model.gtd7d::<f64>(&input);
//! assert!(output.total_mass > 0.0);
//! ```
//!
//! Code that gets its position and time from somewhere with a type never meets this.
//!
//! # Writing against it
//!
//! Three rules, each of which the model's own code obeys throughout:
//!
//! - **Magnitude comparisons branch on [`DualStruct::re`].** `num-dual` does order its
//!   dual numbers, and it does so on the real part alone, so `x > y` would agree with
//!   `x.re() > y.re()`. The model writes the second anyway: every comparison in it — the
//!   altitude bands, the exponent clamps, the guards in `ccor` and `dnet` — is a question
//!   about a real quantity, and saying so at each site is worth the few characters.
//! - **Coefficients stay `f64`.** `T op f64` is available, so lifting three thousand
//!   coefficients to `T` would only cost memory and cache.
//! - **`f64` on the left of an operator does not exist.** `T * 2.0` works; `2.0 * T`
//!   cannot, because of the orphan rule. Write the scalar on the left, or
//!   [`Scalar::lift`] the constant. Where the reference puts a constant on the left of a
//!   *non-commutative* operator, lifting is the only way to keep its arithmetic.
//!
//! # A dual evaluation's value is not bit for bit the `f64` one
//!
//! It agrees to about 1e-14. The difference is `num-dual`'s arithmetic rather than
//! anything this crate does, and two places account for it:
//!
//! - **Division goes through the reciprocal.** `Dual / Dual` is computed as
//!   `self.re * other.re.recip()`, two roundings where a division is one, so a dual
//!   division can land a last bit from the `f64` division. The model divides often — every
//!   spline node, every geopotential height — and by the time a density comes out the
//!   differences have compounded to tens of units in the last place.
//! - **A non-integer power is not the platform `pow`.** `powf(n)` is computed as
//!   `x^(n-3) * x * x * x`, which rounds differently. Helium and argon raise a temperature
//!   ratio to 0.62 and 1.17.
//!
//! Both were measured rather than assumed. Correcting them here was tried and dropped: a
//! real part can be fixed exactly, by adding back a difference that is itself exactly
//! representable, but the same cannot be done for division without replacing every `/` in
//! the model, so the result was a slower model that was still not exact. At `f64` every
//! method is a thin forward to the `std` function of the same name, so that path is the
//! one that matches the reference bit for bit, and it is where the validation rests.

pub use num_dual::{DualNum, DualStruct};
pub use num_traits::Signed;

/// A scalar the model can be evaluated in.
///
/// Satisfied by `f64`, by `num-dual`'s dual numbers over it, and by anything else that
/// implements [`DualNum`] at `f64` precision and is [`Copy`]. There is nothing to
/// implement: the blanket implementation below covers every such type.
///
/// The methods the model computes with are [`DualNum`]'s, [`DualStruct::re`] and
/// [`Signed::abs`]; see the module documentation for what to import.
pub trait Scalar: DualNum<Primitive = f64> + DualStruct<Real = f64> + Copy {
    /// A real constant as a value of this type, with zero derivative.
    ///
    /// `num-dual` spells this `Self::from`, which is a less obvious thing to read in the
    /// middle of an expansion than `T::lift`, and which the model writes a hundred and
    /// twenty times.
    fn lift(x: f64) -> Self {
        Self::from(x)
    }
}

impl<T: DualNum<Primitive = f64> + DualStruct<Real = f64> + Copy> Scalar for T {}

#[cfg(test)]
mod tests {
    use super::*;

    /// The blanket implementation reaches the primitive, the dual types, and a nesting of
    /// them that this crate never names — which is the point of it being blanket.
    #[test]
    fn test_the_blanket_implementation_reaches_what_it_should() {
        fn assert_scalar<T: Scalar>(x: T) -> f64 {
            // One method from each of the four traits, to show they are all available.
            (x.exp().ln().abs() + T::lift(1.0)).re()
        }

        assert_eq!(assert_scalar(2.0_f64), 3.0);
        assert_eq!(assert_scalar(num_dual::Dual64::from(2.0)), 3.0);
        assert_eq!(assert_scalar(num_dual::Dual2_64::from(2.0)), 3.0);
        assert_eq!(assert_scalar(num_dual::HyperDual64::from(2.0)), 3.0);
        // A dual of a dual, which the old list of implementations did not include.
        assert_eq!(
            assert_scalar(num_dual::Dual::<num_dual::Dual64>::from(2.0)),
            3.0
        );
    }

    /// `lift` is a constant: the value with no derivative, whatever the order.
    #[test]
    fn test_lift_carries_no_derivative() {
        use num_dual::Dual64;

        let constant = Dual64::lift(3.5);
        assert_eq!(constant.re, 3.5);
        assert_eq!(constant.eps, 0.0);

        let seeded = Dual64::from(3.5).derivative();
        assert_eq!(seeded.re(), 3.5);
        assert_eq!(seeded.eps, 1.0);
    }

    /// Ordering a dual number compares its real part, which is what the model's branches
    /// want, so the two forms agree.
    #[test]
    fn test_ordering_a_dual_compares_the_real_part() {
        use num_dual::Dual64;

        let small = Dual64::from(1.0).derivative();
        let large = Dual64::from(2.0);
        assert!(small < large);
        assert!(small.re() < large.re());

        // Equal real parts compare equal even though the derivatives differ, where a
        // lexicographic ordering would have made one of them the smaller.
        let other = Dual64::from(1.0);
        assert_ne!(small.eps, other.eps);
        assert_eq!(small.partial_cmp(&other), Some(core::cmp::Ordering::Equal));
    }
}

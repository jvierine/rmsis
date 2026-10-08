#![cfg(feature = "msis21")]
use rmsis::msis21::{Input, Model};

#[test]
#[ignore = "requires original NRL sources/parameters; set RMSIS21_LIBRARY and RMSIS21_PARAMETERS"]
fn linked_backend_matches_pymsis21() {
    let library = std::env::var("RMSIS21_LIBRARY").expect("RMSIS21_LIBRARY");
    let parameters = std::env::var("RMSIS21_PARAMETERS").expect("RMSIS21_PARAMETERS");
    // SAFETY: this test uses the trusted adapter compiled from the supplied sources.
    let model = unsafe { Model::load(&library, &parameters) }.unwrap();
    let second = unsafe { Model::load(&library, &parameters) }.unwrap();
    let cases: serde_json::Value = serde_json::from_str(include_str!("reference21.json")).unwrap();
    // Alternate daily/history modes to force switch reinitialization/cache reset.
    let cases = cases.as_array().unwrap();
    for index in 0..cases.len() {
        let i = if index % 2 == 0 {
            index / 2
        } else {
            cases.len() / 2 + index / 2
        };
        let c = &cases[i];
        let v = |name: &str| c[name].as_f64().unwrap();
        let input = Input {
            day_of_year: v("day_of_year"),
            seconds_of_day: v("seconds_of_day"),
            altitude_km: v("altitude_km"),
            latitude_deg: v("latitude_deg"),
            longitude_deg: v("longitude_deg"),
            f107: v("f107"),
            f107_average: v("f107_average"),
            ap: std::array::from_fn(|j| c["ap"][j].as_f64().unwrap()),
            storm_time: c["mode"] == "history",
        };
        let out = if index % 3 == 0 {
            second.evaluate(&input)
        } else {
            model.evaluate(&input)
        }
        .unwrap();
        for channel in 0..11 {
            let actual = if channel == 10 {
                Some(out.temperature_k)
            } else {
                out.densities[channel]
            };
            let expected = c["expected"][channel].as_f64();
            match (actual, expected) {
                (Some(a), Some(e)) => {
                    let tolerance = if channel == 10 {
                        0.003
                    } else {
                        e.abs() * 5e-5 + if channel == 0 { 1e-35 } else { 1e-3 }
                    };
                    assert!(
                        (a - e).abs() <= tolerance,
                        "case {i} channel {channel}: {a} vs {e}"
                    );
                }
                (None, None) => (),
                _ => panic!("case {i} channel {channel}: missing-value disagreement"),
            }
        }
        assert!(out.exospheric_temperature_k.is_finite());
    }
}

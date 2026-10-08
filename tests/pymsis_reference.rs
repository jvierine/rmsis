use rmsis::{Input, MagneticActivity, Model};

#[test]
fn matches_pymsis_version_zero() {
    let cases: serde_json::Value = serde_json::from_str(include_str!("reference.json")).unwrap();
    let model = Model::new();
    let mut failures = Vec::new();
    for (index, c) in cases.as_array().unwrap().iter().enumerate() {
        let v = |name: &str| c[name].as_f64().unwrap();
        let ap: [f64; 7] = std::array::from_fn(|i| c["ap"][i].as_f64().unwrap());
        let magnetic = if c["mode"] == "daily" {
            MagneticActivity::Daily { ap: ap[0] }
        } else {
            MagneticActivity::History { ap }
        };
        let input = Input::consistent(
            v("day_of_year") as i32,
            v("seconds_of_day"),
            v("altitude_km"),
            v("latitude_deg"),
            v("longitude_deg"),
            v("f107_average"),
            v("f107"),
            magnetic,
        );
        let out = model.gtd7d(&input);
        let actual = [
            out.total_mass,
            out.molecular_nitrogen,
            out.molecular_oxygen,
            out.atomic_oxygen,
            out.helium,
            out.hydrogen,
            out.argon,
            out.atomic_nitrogen,
            out.anomalous_oxygen,
            f64::NAN,
            out.temperature,
        ];
        for (channel, (&a, e)) in actual
            .iter()
            .zip(c["expected"].as_array().unwrap())
            .enumerate()
        {
            if channel == 9 {
                assert!(e.is_null(), "MSISE-00 has no NO");
                continue;
            }
            if let Some(expected) = e.as_f64() {
                // pymsis MSISE-00 evaluates in single precision. Exponentiated
                // upper-atmosphere densities amplify temperature rounding.
                let tolerance = if channel == 10 {
                    0.01
                } else {
                    expected.abs() * 1.5e-4 + if channel == 0 { 1e-35 } else { 1e-3 }
                };
                if !a.is_finite() || (a - expected).abs() > tolerance {
                    failures.push(format!("case {index} channel {channel}: Rust={a} pymsis={expected} relative={} altitude={}",(a-expected)/expected,v("altitude_km")));
                }
            } else {
                assert_eq!(a, 0., "case {index}, undefined species must be zero");
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} discrepancies:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

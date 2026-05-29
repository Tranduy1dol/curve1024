use curve1024::{Curve1024BaseField, Curve1024Config, PrimeFieldConfig, SWCurveConfig};

#[test]
fn test_anomalous_attack_resistance() {
    let p = Curve1024BaseField::MODULUS;
    let r = Curve1024Config::ORDER;
    assert_ne!(r, p);
}

#[test]
fn test_anomalous_weak_curve_simulation() {
    let toy_p: u64 = 11;
    let toy_group_order: u64 = 11;
    let toy_t = toy_p + 1 - toy_group_order;
    assert_eq!(toy_t, 1);
}
